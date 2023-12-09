import os
import sys
import numpy
import pandas
from datetime import datetime
from glob import glob
import metpy.calc as mcalc
from metpy.units import units
from scipy.ndimage import generic_filter
import itertools

from util import get_x_and_y_proj_coords_from_latlon
from grib_reader_ensemble import load_grib,_field_keyval_mapping
from plotmap import PlotMap
import cartopy_util


def get_pydatetime_from_numpy_datetime(dt):
    return pandas.Timestamp(dt).to_pydatetime()

def get_neighborhood_condition(upper_limit,lower_limit):
    def _threshold_any(grid):
        conditional_bool_grid = numpy.where(
            numpy.logical_and(
                grid<=upper_limit,
                grid>=lower_limit
            ),1,0)
        if conditional_bool_grid.any():
            return 1
        else:
            return 0
    return _threshold_any

class Dataset:

    def __init__(self,grib_file_string="",dataset_name="",field_list=list(_field_keyval_mapping.keys()),ensemble_data=True,plotting_domain_fn=None,plot_out_path=None,initialize_plot=True):
  
        self._get_file_list(grib_file_string)
        self.is_ensemble_data = ensemble_data
        self.data = self._read_grib_with_iris(field_list,self.is_ensemble_data)
        self.plot_out_path = self._get_output_path_for_plots(plot_out_path)
        self.cartopy = self._set_cartopy_plotting_metadata(plotting_domain_fn)
        self.dataset_name = dataset_name
        self.res = self.data.res


        self.time_dimension = list(self.data.dims).index('time')
        self.latitude,self.longitude = self._get_coordinate_arrays() 
        
        #get any isobaric DataArray to define isobaric coordinate information
        isobaric_var = [d for d in self.data.data_vars if "pressure" in self.data[d].dims][0]
        self.isobaric_dim_names = self.data[isobaric_var].dims
        self.isobaric_dim_shape = self.data[isobaric_var].shape
        self.isobaric_dim = self.data[isobaric_var].pressure

    
        self.valid_time = [get_pydatetime_from_numpy_datetime(dt) for dt in self.data['time'].values]

        try:
            self.issuance_time = [
                get_pydatetime_from_numpy_datetime(dt) 
                for dt in self.data['forecast_reference_time'].values
            ]
        except TypeError:
            try:
                self.issuance_time = get_pydatetime_from_numpy_datetime(
                    self.data['forecast_reference_time'].values[0]
                )
            except:
                self.issuance_time = None

        if initialize_plot:
            self.plotter = PlotMap(
                self.plot_out_path,
                self.dataset_name,
                self.cartopy,
                self.latitude,
                self.longitude
            )
            
    
    def _set_cartopy_plotting_metadata(self,plotting_domain_fn):
        if plotting_domain_fn:
            cartopy = cartopy_util.PlotCartopy(plotting_domain_fn)
        else:
            cartopy = cartopy_util.PlotCartopy(self.data)
        return cartopy

    def _get_output_path_for_plots(self,plot_out_path):                
        if plot_out_path is None:
            plot_out_path = os.path.join(
                os.path.basename(self._input_grib_files[0]),
                'plots'
            )
        return plot_out_path

    def _get_coordinate_arrays(self):
        latlon_ds = self.data[['latitude','longitude']].reset_coords()
        latlon_ds['longitude'] = latlon_ds.longitude.where(
            latlon_ds.longitude<180.,
            latlon_ds.longitude-360.
        )
        return latlon_ds.latitude.values,latlon_ds.longitude.values

    def _get_file_list(self,nc_file_string):
        if "*" in nc_file_string:
            nc_file_list = glob(nc_file_string)
            nc_file_list.sort()
            if len(nc_file_list) == 0:
                raise FileExistsError("No files found using provided file glob")
            
        elif os.path.isfile(nc_file_string):
            nc_file_list = [nc_file_string]
        else:
            raise FileExistsError(f"{nc_file_string} could not be found")

        self._input_nc_files = nc_file_list

    def _read_nc_with_xarray(self,field_list,ensemble):
        #data = load_grib(self._input_nc_files,field_list=field_list,ensemble=ensemble)
        dataset_list = [xarray.open_dataset(f) for f in self._input_nc_files]
        return dataset_list
    
    def _get_all_dimensions_but_x_and_y(self):
        alldims = self.ds.dims
        drop_dim_dict = {l:0 for l in alldims if l not in ['projection_x_coordinate','projection_y_coordinate']}
        return drop_dim_dict

    def get_profile_dataset(self, inlat:numpy.ndarray, inlon:numpy.ndarray):
        
        points_df = get_x_and_y_proj_coords_from_latlon(inlat,inlon,self.crs)
        
        profile_ds = self.data.sel(
            projection_x_coordinate=points_df['x'].to_xarray(),
            projection_y_coordinate=points_df['y'].to_xarray(),
            method = 'nearest'
        )

        return profile_ds

    def calculate_isobaric_frontogenesis(self,pressure_level):
        potential_temperature = mcalc.potential_temperature(
            numpy.ones_like(self.data['air_temperature_isobaric'].sel(pressure=pressure_level)*pressure_level)*units('hPa'),
            self.data['air_temperature_isobaric'].sel(pressure=pressure_level)*units('K')
        )
        fgen = mcalc.frontogenesis(
            potential_temperature,
            self.data['u_wind_isobaric'].sel(pressure=pressure_level)*units('m/s'),
            self.data['v_wind_isobaric'].sel(pressure=pressure_level)*units('m/s'),
            dx = self.res * units('km'),
            dy = self.res * units('km')
        )
        return fgen

    def append_3D_refl_from_netcdf_files(self,wrf_ensemble_netcdf_path):

        if self.is_ensemble_data:
            looping_bounds = list(itertools.product(self.valid_time,self.data.member.values))
            selectors = [{'time':t,'member':m} for t,m in looping_bounds]
        else:
            selectors = [{'time':t} for t in ds.valid_time]

        self.data['refl_10cm_isobaric'] = (
            self.isobaric_dim_names,
            numpy.ones(self.isobaric_dim_shape)
        )

        for select_dim in selectors:
            valid_time = select_dim.get('time')
            member = select_dim.get('member')
            wrf_fn = os.path.join(
                wrf_ensemble_netcdf_path,
                f"wrfout_d01_{valid_time:%Y-%m-%d_%H:%M:%S}"
            )
            if member:
                wrf_fn += f"_{member:04}"
            
            wrf_ds = xarray.open_dataset(wrf_fn)
            pressure = (wrf_ds.P + wrf_ds.PB).squeeze()/100.
            refl_10cm = wrf_ds.REFL_10CM.squeeze() * metpy.units.units('')
            
            isobaric_refl = metpy.interpolate.log_interpolate_1d(
                self.data.pressure.to_numpy(),
                pressure,
                refl_10cm,
                axis=0
            )
            self.data['refl_10cm_isobaric'].loc[select_dim] = isobaric_refl
            
    def calculate_NMEP(self,field_name,lower_threshold,upper_threshold):
        field_ds = self.data[field_name].copy(deep=True)
        for member in field_ds.member:
            field_ds.loc[dict(member=member)] = generic_filter(
                field_ds.sel(member=member).to_numpy(),
                get_neighborhood_condition(
                    upper_threshold,
                    lower_threshold
                ),
                size = 4,
                mode = 'nearest'
            )
        
        nmep_ds = field_ds.mean(dim='member') * 100
        return nmep_ds

    def calculate_PMM(self,field_name):
        field_ds = self.data[field_name].copy(deep=True)

        assert 'member' in field_ds.dims
        field_ds_mean = field_ds.mean(dim='member')
        field_ds_mean.name = field_ds_mean.name + "_pmm"
        nmems = field_ds.member.shape[0]

        for time in field_ds.time:
            ensemble_mean = field_ds_mean.sel(time=time)
            ensemble_mean_value_indices = numpy.argsort(ensemble_mean.to_numpy().flatten())
            ensemble_full_value_magnitudes = numpy.sort(field_ds.sel(time=time).to_numpy().flatten())[0::nmems]
            
            pmm_flat = numpy.ones_like(ensemble_mean_value_indices)
            pmm_flat[ensemble_mean_value_indices] = ensemble_full_value_magnitudes
            pmm = pmm_flat.reshape(ensemble_mean.shape)
            
            field_ds_mean.loc[dict(time=time)] = pmm

        return field_ds_mean 
                 
    
    def plot_contourfilled_map(self,field_name,selector_dict,ax=None,**kwargs):
        prop = kwargs.pop('prop',{})
        map_prop = kwargs.pop('map_prop',{})
        
        data = self.data.sel(selector_dict)[field_name].copy(deep=True)
        
        self.plotter.init_plot()
        self.plotter.set_plot_time_title(
            pandas.to_datetime(data.time.values.item()),
            issuance_time=pandas.to_datetime(data.forecast_reference_time.values.item())
        )
        self.plotter._set_plot_custom_title(figure_title=f'{field_name} [color-filled]')
    
        self.plotter.define_colormap()
        self.plotter._configure_colormap_legend(
            legend_name=f"{field_name}      ",
            bounds=numpy.arange(data.min(),data.max(),data.min()),
            ticks=numpy.arange(data.min(),data.max(),data.min())
        )
        contourfilled = self.plotter.add_filled_contour(
            data.to_numpy(),
            numpy.arange(data.nanmin(),data.nanmax,data.nanmin()),
            self.plotter.cmap
        )
        self.plotter._set_xylims()
        self.plotter._set_plot_filename(f"{field_name}",data.time.values.item(),issuance_time=data.forecast_reference_time.values.item())
        self.plotter.save_plot()

    def plot_surface_field(self,field_name,selector_dict={},ax=None,cartopy_proj=None,save_path=None,**kwargs):


        #Retrieve kwargs
        prop = kwargs.pop('prop',{})
        map_prop = kwargs.pop('map_prop',{})

        try:
            self.plot_obj
        except:
            self.plot_obj = PlotMap()

        if cartopy_proj:
            self.plot_obj.proj = cartopy_proj
        else:
            self.plot_obj.create_cartopy(proj='PlateCarree',central_longitude=0.0)

        data = self.surface_data.get(field_name).sel(selector_dict)
        
        #if data == None:
        #    print(f'{field_name} cannot be found in surface dataset. Please investigate further')
        
        plot_ax = self.plot_obj.plot_map(data,self.latitude,self.longitude,ax=ax,save_path=save_path)
        self.plot_obj.ax.contourf(self.longitude,self.latitude,data.values)

    #def plot_MSLP_GH500mb_CREFL(self,plot_prefix_name=None,subset_domain_grid_file=None):

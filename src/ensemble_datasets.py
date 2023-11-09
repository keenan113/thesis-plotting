import os
import sys
import numpy
import pandas
from datetime import datetime
from glob import glob
import metpy.calc as mcalc
from metpy.units import units

from util import get_x_and_y_proj_coords_from_latlon
from grib_reader_ensemble import load_grib,_field_keyval_mapping
from plotmap import PlotMap


def get_pydatetime_from_numpy_datetime(dt):
    return pandas.Timestamp(dt).to_pydatetime()

class Dataset:

    def __init__(self,grib_file_string="",field_list=list(_field_keyval_mapping.keys()),ensemble_data=True):
  
        self._get_file_list(grib_file_string)

        self._read_grib_with_iris(field_list,ensemble_data)

        self.proj = None

    def _get_file_list(self,grib_file_string):
        if "*" in grib_file_string:
            grib_file_list = glob(grib_file_string)
            grib_file_list.sort()
            if len(grib_file_list) == 0:
                raise FileExistsError("No files found using provided file glob")
            
        elif os.path.isfile(grib_file_string):
            grib_file_list = [grib_file_string]
        else:
            raise FileExistsError(f"{grib_file_string} could not be found")

        self._input_grib_files = grib_file_list

    def _read_grib_with_iris(self,field_list,ensemble):
        grib_datasets = load_grib(self._input_grib_files,field_list=field_list,ensemble=ensemble)
        
        try:
            assert 'surface-dataset' in grib_datasets.keys()
            self.surface_data = grib_datasets['surface-dataset']
            varname = list(self.surface_data.data_vars)[0]
        except AssertionError:
            raise ValueError("Expected at least the surface dataset...not found after load")
        
        try:
            self.isobaric_data = grib_datasets['isobaric-dataset']
        except KeyError:
            self.isobaric_data = None

        try:
            self.hybrid_data = grib_datasets['hybrid-dataset']
        except KeyError:
            self.hybrid_data = None

        self.crs = self.surface_data[varname].crs

        alldims = self.surface_data[varname].dims
        drop_dim_dict = {l:0 for l in alldims if l not in ['projection_x_coordinate','projection_y_coordinate']}
        
        
        #if 'time' in self.surface_data.dims:
        self.time_dimension = list(self.surface_data.dims).index('time')
        self.latitude = self.surface_data.isel(drop_dim_dict)['latitude'].values
        self.longitude = self.surface_data.isel(drop_dim_dict)['longitude'].where(
            self.surface_data.isel(drop_dim_dict)['longitude']<180,
            self.surface_data.isel(drop_dim_dict)['longitude']-360.
        ).values
        self.valid_time = [get_pydatetime_from_numpy_datetime(dt) for dt in self.surface_data['time'].values]
        self.issuance_time = get_pydatetime_from_numpy_datetime(self.surface_data['forecast_reference_time'].values[0])
        #else:
        #    self.time_dimension = None
        #    self.latitude = self.surface_data['latitude'].isevalues
        #    self.longitude = self.surface_data['longitude'].where(
        #        self.surface_data['longitude']<180,
        #        self.surface_data['longitude']-360.
        #    ).values
        #    self.valid_time = get_pydatetime_from_numpy_datetime(self.surface_data['time'].values)
        #    self.issuance_time = get_pydatetime_from_numpy_datetime(self.surface_data['forecast_reference_time'].values)

    def get_isobaric_profile_dataset(self, inlat:numpy.ndarray, inlon:numpy.ndarray):
        if not self.isobaric_data:
            raise ValueError("No isobaric data loaded - so outputting a profile is impossible")
        
        points_df = get_x_and_y_proj_coords_from_latlon(inlat,inlon,self.crs)
        
        profile_ds = self.isobaric_data.sel(
            projection_x_coordinate=points_df['x'].to_xarray(),
            projection_y_coordinate=points_df['y'].to_xarray(),
            method = 'nearest'
        )

        return profile_ds

    def calculate_isobaric_frontogenesis(self,pressure_level):
        potential_temperature = mcalc.potential_temperature(
            numpy.ones_like(self.isobaric_data['temperature'].sel(level=pressure_level)*pressure_level)*units('hPa'),
            self.isobaric_data['temperature'].sel(level=pressure_level)*units('Kelvin')
        )
        fgen = mcalc.frontogenesis(
            potential_temperature,
            self.isobaric_data['u'].sel(level=pressure_level)*units('m/s'),
            self.isobaric_data['v'].sel(level=pressure_level)*units('m/s'),
            #dx=
        )

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





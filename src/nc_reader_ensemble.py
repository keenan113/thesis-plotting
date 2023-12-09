import os
import warnings
import numpy
import iris
import pandas
import xarray
import cartopy.crs as ccrs
import xarray
import numpy
import util
import cartopy_util


def _get_member_from_filename(f):
    if len(os.path.basename(f).split('_'))==5:
        member = int(os.path.basename(f).split('.')[0].split('_')[-1])
    else:
        member = None
    return member

def _open_single_file_dataset_ensemble(fn):
    ds = xarray.open_dataset(fn)
    ds = ds.swap_dims({"Time":"XTIME"})
    
    member = _get_member_from_filename(fn)
    if member: 
        ds = ds.assign_coords({"member":member}).expand_dims("member")

    return ds    

def load_nc_files(filenames,ensemble=False):

    data_array_list = [_open_single_file_dataset_ensemble(fn) for fn in filenames]

    dataset = xarray.combine_by_coords(data_array_list,combine_attrs='drop_conflicts')

    dataset.attrs = {"crs":cartopy_util.PlotCartopy(filenames[0])}

    return dataset


def _generate_projected_coordinate_data_arrays(da):
    pc = ccrs.PlateCarree()
    ys,xs = numpy.meshgrid(da.projection_y_coordinate,da.projection_x_coordinate)
    longitude,latitude,_ = pc.transform_points(da.crs,xs,ys).T
    longitude_da = _create_latlon_data_arrays(da,longitude)
    longitude_da.name = 'longitude'
    latitude_da = _create_latlon_data_arrays(da,latitude)
    latitude_da.name = 'latitude'
    return longitude_da,latitude_da


def _create_latlon_data_arrays(da,ll_array):
    coords = {
        'projection_y_coordinate':da['projection_y_coordinate'],
        'projection_x_coordinate':da['projection_x_coordinate']
    }
    dims = ['projection_y_coordinate','projection_x_coordinate']
    return xarray.DataArray(data=ll_array,coords=coords,dims=dims)

def _rename_variables(inds):
    """
    Renames fields in a given dataset to a single set of standard names
    that make matching datasets simple.

    Args:
        inds (xarray.Dataset): dataset containing any number of meteorological fields
    Returns:
        (xarray.Dataset): dataset with fields renamed to standard names
        
        'air_temperature_surface':'temperature',
        'dewpoint_temperature_surface':'dewpoint',
        'u_wind_surface':'uwind',
        'v_wind_surface':'vwind',
        'pressure_surface':'pressure',
        'terrain_height_surface':'height',
        'specific_humidity_isobaric':'specific_humidity',
        'air_temperature_isobaric':'temperature',
        'u_wind_isobaric':'uwind',
        'v_wind_isobaric':'vwind',
        'geopotential_height_isobaric':'height'
    """

    rename_dict = {
        "mslp_MAPS_surface":"mslp_surface"
    }
    dataset_variables = list(inds.data_vars)
    rename_dict = {key:value for key,value in rename_dict.items() if key in dataset_variables}
    return inds.rename(rename_dict)

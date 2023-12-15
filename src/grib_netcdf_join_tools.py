import numpy
import xarray
import metpy

def add_netcdf_field_to_grib_dataset_isobaric(grib_ds,netcdf_ds,netcdf_field_name):
    assert (grib_ds.member == netcdf_ds.member).all() and (grib_ds.time.values == netcdf_ds.XTIME.values).all()
    new_field_isobaric = metpy.interpolate.log_interpolate_1d(
        grib_ds.pressure.to_numpy(),
        (netcdf_ds.P + netcdf_ds.PB)/100.,
        netcdf_ds[netcdf_field_name] * metpy.units.units(""),
        axis=netcdf_ds.P.dims.index('bottom_top')
    )
    new_field_name = netcdf_ds[netcdf_field_name].name.lower() + "_isobaric" 
    netcdf_ds = netcdf_ds.rename({'XTIME':'time','bottom_top':'pressure',
                        'south_north':'projection_y_coordinate',
                        'west_east':'projection_x_coordinate'})
    grib_ds[new_field_name] = ( netcdf_ds[netcdf_field_name].dims,
                                new_field_isobaric.magnitude )
    return grib_ds

def _get_netcdf_file_list(inds):
    times = inds.data.time
    if 'member' in inds.data.dims:
        members = inds.data.member
    
    return   
    
    

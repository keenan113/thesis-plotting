import numpy
import xarray
import cartopy

def _define_projection_extents(dest_crs,geobounds):
    # Need to modify the extents for the new projection
    
    src_crs = cartopy.crs.PlateCarree()
    xs, ys, _ = dest_crs.transform_points(
        src_crs,
        numpy.array([geobounds['bottom_left']['lon'], geobounds['top_right']['lon']]),
        numpy.array([geobounds['bottom_left']['lat'], geobounds['top_right']['lat']])
    ).T
    _xlimits = xs.tolist()
    _ylimits = ys.tolist()

    return (_xlimits, _ylimits)

def _get_cartopy_crs(wrfds):
    WRF_EARTH_RADIUS = 6370000.

    globe = cartopy.crs.Globe(
        ellipse=None,
        semimajor_axis=WRF_EARTH_RADIUS,
        semiminor_axis=WRF_EARTH_RADIUS,
        nadgrids="@null"
    )

    cutoff = -30.0 if wrfds.attrs['MOAD_CEN_LAT'] >= 0 else 30.0

    ccrs = cartopy.crs.LambertConformal(
        central_longitude = wrfds.attrs['STAND_LON'],
        central_latitude = wrfds.attrs['MOAD_CEN_LAT'],
        standard_parallels = (wrfds.attrs['TRUELAT1'],wrfds.attrs['TRUELAT2']),
        globe = globe,
        cutoff = cutoff
    )
    return ccrs

def _define_grid_bounds(lats,lons):
    #TODO replace this simple function with fully fleshed out data classes
    geobounds = {
        "bottom_left":{
            "lon":lons[0,0],
            "lat":lats[0,0]
        },
        "top_right":{
            "lon":lons[-1,-1],
            "lat":lats[-1,-1]
        }, 
    }
    return geobounds

def get_cartopy_metadata(sample_netcdf_grid_fn):

    ds = xarray.open_dataset(sample_netcdf_grid_fn)
    lats = ds.XLAT.values[0,:,:]
    lons = ds.XLONG.values[0,:,:]
    geobounds = _define_grid_bounds(lats,lons)
    projection = _get_cartopy_crs(ds)
    xlim, ylim = _define_projection_extents(projection,geobounds)
    cartopy_metadata = {
        "lats":lats,
        "lons":lons,  
        "geobounds":geobounds,
        "projection":projection,
        "xlim":xlim,
        "ylim":ylim
    }
    return cartopy_metadata

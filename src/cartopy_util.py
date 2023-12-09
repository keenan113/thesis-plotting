import numpy
import xarray
import cartopy
from geobound_util import Geobound,CoordPair

    
def _define_projection_extents(dest_crs,geobounds):
    # Need to modify the extents for the new projection
    
    src_crs = cartopy.crs.PlateCarree()
    xs, ys, _ = dest_crs.transform_points(
        src_crs,
        numpy.array([geobounds.bottom_left.longitude, geobounds.top_right.longitude]),
        numpy.array([geobounds.bottom_left.latitude, geobounds.top_right.latitude])
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

def get_cartopy_metadata_from_wrfnetcdf(sample_netcdf_grid_fn):

    ds = xarray.open_dataset(sample_netcdf_grid_fn)
    lats = ds.XLAT.values[0,:,:]
    lons = ds.XLONG.values[0,:,:]
    projection = _get_cartopy_crs(ds)
    geobounds = Geobound(latitudes=lats,longitudes=lons,crs=projection) 
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

class PlotCartopy:
    def __init__(self,inds):
        if type(inds) is str: 
            cartopy_meta = get_cartopy_metadata_from_wrfnetcdf(inds)
            self.projection = cartopy_meta['projection']
            self.geobounds = cartopy_meta['geobounds']
            self.xlim = cartopy_meta['xlim']
            self.ylim = cartopy_meta['ylim']
            self.xlats = cartopy_meta['lats']
            self.xlons = cartopy_meta['lons'] 

        elif type(inds) is xarray.Dataset:
            self.projection = inds.crs
            self.geobounds = Geobound(latitudes=inds.latitude,longitudes=inds.longitude,crs=self.projection)
            self.xlim, self.ylim = _define_projection_extents(self.projection,self.geobounds)
        else:
            print('Unknown input file format (either filename or ensemble dataset should be provided)')
    
    def _get_xy_plot_limits(self):
        if self.geobounds.bottom_left.x is not None and self.geobounds.top_right.x is not None:
            self.xlim = [self.geobounds.bottom_left.x, self.geobounds.top_right.x]
        else:
            self.xlim = None
        if self.geobounds.bottom_left.y is not None and self.geobounds.top_right.y is not None:
            self.ylim = [self.geobounds.bottom_left.y, self.geobounds.top_right.y]
        else:
            self.ylim = None

        if self.xlim is None or self.ylim is None:
            self.xlim, self.ylim = _define_projcetion_extents(self.projection, self.geobounds)

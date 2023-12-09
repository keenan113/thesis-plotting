import util
import numpy
import cartopy.crs as ccrs

class Geobound:
    def __init__(self,bottom_left=None,top_right=None,latitudes=None,longitudes=None,crs=None):
        if bottom_left is not None and top_right is not None:
            self.bottom_left = bottom_left
            self.top_right = top_right
            if bottom_left.latitude is None or bottom_left.longitude is None:
                raise ValueError("'bottom_left' does not contain a latitude or longitude attribute")
            if top_right.latitude is None or top_right.longitude is None:
                raise ValueError("'top_right' does not contain a latitude or longitude attribute")
        elif latitudes is not None and longitudes is not None:
            if latitudes.ndim != 2 or longitudes.ndim != 2:
                raise ValueError("'latitude' and 'longitude' arrays must be 2-dimensional."
                                 " Got {latitudes.ndim} dims for 'latitude' and"
                                 " {longitudes.ndim} dims for 'longitude'")
            self.bottom_left = CoordPair(latitude=latitudes[0,0],longitude=longitudes[0,0],crs=crs)
            self.top_right = CoordPair(latitude=latitudes[-1,-1],longitude=longitudes[-1,-1],crs=crs)
        else:
            raise ValueError("Must specify either 'bottom_left' & 'top_right' parameters"
                             " or specify 'latitude' and 'longitude' array parameters")

class CoordPair:
    def __init__(self,latitude=None,longitude=None,x=None,y=None,crs=None):
        
        self.latitude = latitude
        self.longitude = longitude
        self.x = x
        self.y = y
        self.crs = crs

        if self.x is None and self.y is None:
            self._find_projected_xy_values()

        if self.latitude is None and self.longitude is None:
            self._find_latlon()
        
    def _find_projected_xy_values(self): 
        try:
            x,y = util.reproject_latlon_to_xy(self.latitude,self.longitude,self.crs)
            self.x = x
            self.y = y
        except:
            self.x = None
            self.y = None

    def _find_latlon(self):
        try:
            lat,lon = util.reproject_xy_to_latlon(self.x,self.y,self.crs)
            self.latitude = lat
            self.longitude = lon
        except:
            self.latitude = None
            self.longitude = None

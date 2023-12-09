import pandas
import numpy
import cartopy.crs as ccrs

def get_x_and_y_proj_coords_from_latlon(latitude,longitude,data_crs):
    """
    Args:
        latitude (numpy.ndarray)
        longitude (numpy.ndarray)
        data_crs (cartopy.CRS)
    Returns:
        coordinate_dataframe (pandas.DataFrame)
    """
    out_coordinates = data_crs.transform_points(
        ccrs.PlateCarree(),
        x=longitude,
        y=latitude
        )

    coord_dataframe = pandas.DataFrame({
        "x":[c[0] for c in out_coordinates],
        "y":[c[1] for c in out_coordinates]
    })

    return coord_dataframe

def reproject_latlon_to_xy(latitude,longitude,src_crs,dest_crs = ccrs.PlateCarree()):

    if type(latitude) is not numpy.ndarray or type(longitude) is not numpy.ndarray:
        latitude = numpy.array(latitude)
        longitude = numpy.array(longitude)

    x,y,_ = src_crs.transform_points(
        dest_crs,
        x = longitude,
        y = latitude
    ).T
    return x,y

def reproject_xy_to_latlon(x,y,dest_crs):

    if type(x) is not numpy.ndarray or type(y) is not numpy.ndarray:
        x = numpy.array(x)
        y = numpy.array(y)

    latitude,longitude,_ = ccrs.PlateCarree().transform_points(
        dest_crs,
        x = x,
        y = y
    ).T
    return latitude,longitude

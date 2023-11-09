import pandas
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
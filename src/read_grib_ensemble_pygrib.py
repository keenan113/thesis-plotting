import pygrib
import numpy as np
import cartopy.crs as crs
from iris_grib._load_convert import unscale, ellipsoid_geometry, ellipsoid
import iris.coord_systems as icoord_systems

def get_metadata(grbmsg):
    metadict = {
    'issuance_date':grbmsg.analDate,
    'valid_date':grbmsg.validDate,
    }
    _,xlats,xlons = grbmsg.data()
    projection = _get_grid_coordinate_reference_system(grbmsg)#crs.PlateCarree()
    metadict['xlats'] = xlats
    metadict['xlons'] = xlons
    metadict['projection'] = projection
    return metadict

def read_in_ensemble_variable(list_of_files,keyword_args):
    """
    This function will read a list of files, assuming that each file represents
    a grib file with wrfout fields at a single time. Each entry in the list is
    then a different member in an ensemble. In addition, this function will only
    in a single field from each file.

    Args:
        list_of_files (list): list of absolute paths to wrfout grib2 data with
            each entry representing a different member of an ensemble
        keyword_args (dict): a dictionary of grib2 key/values for obtaining a single
            field from each grib file
    Returns:
        (dictionary): 
            data: a dictionary with the data as a numpy array
            meta: a dictionary that contains a collection of grid data necessary for plotting
    """
    list_of_files.sort()
    grblist = []
    for fn in list_of_files:
        grb = pygrib.open(fn)
        arr = grb.select(**keyword_args)[0]
        grblist.append(arr)
        grb.close()
    meta = get_metadata(grblist[0])
    datalist = [g.values for g in grblist]
    data = np.zeros([datalist[0].shape[0],datalist[0].shape[1],len(datalist)])
    for mem in range(data.shape[2]):
        data[:,:,mem] = grblist[mem].values
    return {'data':data,'meta':meta}

def _get_grid_coordinate_reference_system(section):

    _GRID_ACCURACY_IN_DEGREES = 1e-6

    major, minor, radius = ellipsoid_geometry(section)
    geog_cs = ellipsoid(section['shapeOfTheEarth'], major, minor, radius)
    central_latitude = section['LaD'] * _GRID_ACCURACY_IN_DEGREES
    central_longitude = section['LoV'] * _GRID_ACCURACY_IN_DEGREES
    false_easting = 0
    false_northing = 0
    secant_latitudes = (section['Latin1'] * _GRID_ACCURACY_IN_DEGREES,
                        section['Latin2'] * _GRID_ACCURACY_IN_DEGREES)

    cs = icoord_systems.LambertConformal(central_latitude,
                                        central_longitude,
                                        false_easting,
                                        false_northing,
                                        secant_latitudes=secant_latitudes,
                                        ellipsoid=geog_cs).as_cartopy_crs()
    

    return cs
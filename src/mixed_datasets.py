import os
import numpy
import pandas
from datetime import datetime
from glob import glob

from src import grib_reader,read_grib_ensemble_pygrib

#DEPRECATED

def get_pydatetime_from_numpy_datetime(dt):
    return pandas.Timestamp(dt).to_pydatetime()

class Dataset:

    def __init__(self,grib_file_string="",data_type="iris",variable_dict=None):

        if "*" in grib_file_string:
            grib_file_list = glob(grib_file_string)
            grib_file_list.sort()
            if len(grib_file_list) == 0:
                raise FileExistsError("No files found using provided file glob")
        elif os.path.isfile(grib_file_string):
            grib_file_list = [grib_file_string]
        else:
            raise FileExistsError(f"{grib_file_string} could not be found")

        if data_type not in ["iris","pygrib"]:
            raise ValueError("data_type parameter must be either iris or pygrib")
        
        if data_type=='iris':
            self._read_grib_with_iris(grib_file_list,variable_dict)
        else:
            self._read_grib_with_pygrib(grib_file_list,variable_dict)
    
    def _read_grib_with_iris(self,gribfn_list,keyvals):
        datadict = grib_reader.load_grib(gribfn_list)
        assert set(keyvals.keys()).issubset(set(['dataset','varname']))
        ds = datadict[keyvals['dataset']]
        if keyvals['dataset'] == 'isobaric-dataset':
            assert 'level' in keyvals.keys()
            ds = ds.sel(pressure=keyvals['level'])
        
        self.data = ds[keyvals['varname']]
        self.proj = ds[keyvals['varname']].crs

        if 'time' in self.data.dims:
            self.time_dimension = list(self.data.dims).index('time')
            self.latitude = ds['latitude'].isel(time=0).values
            self.longitude = ds['longitude'].isel(time=0).values
            self.valid_time = [get_pydatetime_from_numpy_datetime(dt) for dt in ds['time'].values]
            self.issuance_time = get_pydatetime_from_numpy_datetime(ds['forecast_reference_time'])
        else:
            self.time_dimension = None
            self.latitude = ds['latitude'].values
            self.longitude = ds['longitude'].values
            self.valid_time = get_pydatetime_from_numpy_datetime(ds['time'].values)
            self.issuance_time = get_pydatetime_from_numpy_datetime(ds['forecast_reference_time'].values)            

    def _read_grib_with_pygrib(self,gribfn_list,keyvals):
        datadict = read_grib_ensemble_pygrib.read_in_ensemble_variable(
            gribfn_list,
            keyvals
        )
        self.data = datadict['data']
        self.proj = datadict['meta']['projection']
        self.latitude = datadict['meta']['xlats']
        self.longitude = datadict['meta']['xlons']
        self.valid_time = datadict['meta']['valid_date']
        self.issuance_time = datadict['meta']['issuance_date']
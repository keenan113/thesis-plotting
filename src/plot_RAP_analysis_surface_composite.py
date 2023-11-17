#!/usr/localbin/python3
import argparse
import os
import sys
import pandas
import numpy as np
import itertools
import cartopy_util 
import grib_reader
import plotmap

def main():
    args = parse_args()
    
    cartopy_metadata = cartopy_util.get_cartopy_metadata(args.sample_grid_file)
    
    ds = grib_reader.load_grib(args.infile)    
    latitudes = ds.latitude.values
    longitudes = ds.longitude.values
    wx_plotter = plotmap.PlotMap(
        args.outdir,
        "RAP_analysis",
        cartopy_metadata,
        latitudes,
        longitudes
    )
    
    for time,issuance_time in itertools.product(ds.time,ds.forecast_reference_time):

        MSLP = ds.mslp_MAPS_surface.sel(
            time=time,
            forecast_reference_time=issuance_time
        ).values/100

        CREFL = ds.composite_reflectivity_surface.sel(
            time=time,
            forecast_reference_time=issuance_time
        ).values

        GH500mb = ds.geopotential_height_isobaric.sel(
            time=time,
            forecast_reference_time=issuance_time,
            pressure=500.
        ).values

        valid_time = pandas.to_datetime(time.values)
        
        wx_plotter.plot_MSLP_GH500mb_CREFL(
            valid_time,
            MSLP,
            GH500mb,
            CREFL
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description = "plot RAP analysis"
    )
    parser.add_argument(
        '--infile',
        required = True,
        help = 'Path to a single file for a file glob'
    )
    parser.add_argument(
        '--outdir',
        default = os.path.join(os.getcwd(),'rap_surface_analyses_composites'),
        help = 'path where to write the plots'
    )
    parser.add_argument(
        '--sample-grid-file',
        #default = '/F0/keenan/ThesisProject/Feb2013Case/YesradRun/201302081815/assim_20130208_234500/wrfout_d01_2013-02-08_23:45:00_0001',
        default = '/F0/keenan/ThesisProject/Feb2013Case/2013020406/assim_2013020806/wrfout_d01_2013-02-08_00:00:00_0001',
        help = 'path to any WRF file that contains to projection information that we want to use to define the domain for plotting'
    )
    return parser.parse_args()


if __name__ == '__main__':
    main()

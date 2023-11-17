#!/usr/localbin/python3
import argparse
import os
import sys
import pandas
import numpy as np
import cartopy_util 
import grib_reader
import plotmap

def main():
    args = parse_args()
    
    cartopy_metadata = cartopy_util.get_cartopy_metadata(args.sample_grid_file)
    
    ds = grib_reader.load_grib(args.infile)    
    latitudes = ds['surface-dataset'].latitude.values
    longitudes = ds['surface-dataset'].longitude.values
    wx_plotter = plotmap.PlotMap(
        args.outdir,
        "RAP_analysis",
        cartopy_metadata,
        latitudes,
        longitudes
    )
    
    times = ds['isobaric-dataset'].time
    for time in times.values:

        MSLP, CREFL, GH500mb = preprocess_fields(ds,time)
        valid_time = pandas.to_datetime(time)
        
        wx_plotter.plot_MSLP_GH500mb_CREFL(
            valid_time,
            MSLP,
            GH500mb,
            CREFL
        )

def preprocess_fields(ds,time):
    isobaric_ds = ds['isobaric-dataset']
    surface_ds = ds['surface-dataset']

    MSLP = surface_ds.mslp_MAPS_surface.sel(time=time).values/100
    CREFL = surface_ds.composite_reflectivity_surface.sel(time=time).values
    GH500mb = isobaric_ds.height.sel(time=time,pressure=500.).values
    return MSLP,CREFL,GH500mb

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

#!/usr/localbin/python3
import argparse
import os
import sys
import pandas
import numpy as np
import itertools
from datetime import datetime
import ensemble_datasets
import plotmap

def main():
    args = parse_args()

    os.makedirs(args.outdir,exist_ok=True) 

    ds = ensemble_datasets.Dataset(
        args.infile,
        dataset_name=args.dataset_name,
        ensemble_data=args.ensemble,
        plot_out_path=args.outdir,
        plotting_domain_fn=args.sample_grid_file
    )
    print(f'{datetime.now():%Y%m%dT%H%M}.......Loaded data to dataset')
    if args.pmm:
        crefl_pmm = ds.calculate_PMM('composite_reflectivity_surface')
        ds.data = ds.data.mean(dim='member').assign({'composite_reflectivity_surface':crefl_pmm})
    else:
        ds.data = ds.data.mean(dim='member')

    print(f'{datetime.now():%Y%m%dT%H%M}.......Calculated ensemble mean')

    for valid_time in ds.valid_time:
        MSLP = ds.data.mslp_surface.sel(
            time=valid_time,
        ).values/100
        CREFL = ds.data.composite_reflectivity_surface.sel(
            time=valid_time,
        ).values

        GH500mb = ds.data.geopotential_height_isobaric.sel(
            time=valid_time,
            pressure=500.
        ).values

        
        ds.plotter.plot_MSLP_GH500mb_CREFL(
            valid_time,
            MSLP,
            GH500mb,
            CREFL
        )

        print(f'{datetime.now():%Y%m%dT%H%M}.......Saved plot for {valid_time:%Y%m%dT%H%M}')

def parse_args():
    parser = argparse.ArgumentParser(
        description = "Plot MSLP, 500mb GH, and CREFL"
    )
    parser.add_argument(
        '--dataset-name',
        help = 'name to give dataset (used in plots etc.)',
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
        #default = '/F0/keenan/ThesisProject/Feb2013Case/2013020406/assim_2013020806/wrfout_d01_2013-02-08_00:00:00_0001',
        default = None,
        help = 'path to any WRF file that contains to projection information that we want to use to define the domain for plotting'
    )
    parser.add_argument(
        '--ensemble',
        help = 'flag to indicate data is an ensemble',
        action='store_true'
    )
    parser.add_argument(
        '--pmm',
        action = 'store_true',
        help = 'flag to indicate if a PMM should be used for the contour-filled field. Otherwise a traditional mean will be used.'
    )
    return parser.parse_args()


if __name__ == '__main__':
    main()

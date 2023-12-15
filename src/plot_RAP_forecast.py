#!/usr/localbin/python3

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from matplotlib.cm import get_cmap
import pygrib
import numpy as np
from datetime import datetime, timedelta
import sys
import os
import metpy.calc as mcalc
from metpy.units import units
import pint
from glob import glob
from metpy.plots import ctables
import argparse
from scipy.ndimage import generic_filter

import numpy
import os
import pandas
from src.grib_reader import load_grib



def main():
    args = parse_args()
    os.makedirs(args.output_directory,exist_ok=True)
    rap_grb = args.
    
    
    print(f'Working on time: {data_dict["meta"]["valid_date"]:%Y%m%dT%H%M%S}')

    for mem in range(data_dict['data'].shape[2]):
        data_dict['data'][:,:,mem] = generic_filter(
            data_dict['data'][:,:,mem],
            get_neighborhood_condition(
                args.upper_threshold,
                args.lower_threshold
            ),
            size = 4,
            mode = 'nearest'
        )
            
    data_dict['data'] = data_dict['data'].mean(axis=2)*100
    figname = os.path.join(
        args.output_directory,
        f'WRFEnsemble_NMEP_1kmRefl_between'
        f'{args.lower_threshold}-'
        f'{args.upper_threshold}dBZ'
        f'_{args.name}'
    )
    title = f'WRF-{args.name} Ensemble NMEP 1km AGL Reflectivity | probability between {args.lower_threshold}dBZ and {args.upper_threshold}dBZ in 12km neighborhood'

    plot_probability_reflectivity_field(
        figname,
        title,
        data_dict
    )

def get_neighborhood_condition(upper_limit,lower_limit):
    def _threshold_any(grid):
        conditional_bool_grid = np.where(
            np.logical_and(
                grid<=upper_limit,
                grid>=lower_limit
            ),1,0)
        if conditional_bool_grid.any():
            return 1
        else:
            return 0
    return _threshold_any
    
def get_metadata(grbmsg):
    metadict = {
    'issuance_date':grbmsg.analDate,
    'valid_date':grbmsg.validDate,
    }
    _,xlats,xlons = grbmsg.data()
    projection = crs.PlateCarree()
    metadict['xlats'] = xlats
    metadict['xlons'] = xlons
    metadict['projection'] = projection
    return metadict

def read_in_ensemble_variable(list_of_files,keyword_args):
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
        
def plot_probability_reflectivity_field(figurename, title, datadict):

    fmtstr = '%3.0f'

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(2, 1, 1, projection=datadict['meta']['projection'], anchor='S')
    ax1.set_position([0.05, 0.09, 0.9, 0.9], which='both')
    #ax1.outline_patch.set_linewidth(1.0)
    #ax1.outline_patch.set_zorder(6)
    xlons = datadict['meta']['xlons']
    xlats = datadict['meta']['xlats']
    ax1.set_extent([xlons[0,0],xlons[-1,-1],xlats[0,0],xlats[-1,-1]],datadict['meta']['projection'])
    states = NaturalEarthFeature(
        category='cultural',
        scale='10m',
        facecolor='none', 
        name='admin_1_states_provinces_shp'
    )
    lakes = NaturalEarthFeature(
        category='physical',
        scale='10m',
        facecolor='none',
        name='lakes',
        edgecolor='#ffffff'
    )
    ax1.add_feature(states, edgecolor='#404040',
                    linewidth=1.0, facecolor='none',
                    zorder=4)
    ax1.add_feature(lakes, edgecolor='#404040',
                    linewidth=1.0, zorder=4)
    ax1.set_title(f"Valid: {datadict['meta']['valid_date']:%m/%d/%Y @ %H:%M} UTC",
                    fontsize='large',loc='right')
    ax1.set_title(title, fontsize='medium', loc='left')

    #Import the NWS reflectivity colormap
    cmap = get_cmap('viridis')
    #cmap = ctables.registry.get_colortable('NWSReflectivity')
    #cmap.set_under('white')

    ax2 = fig.add_subplot(2, 1, 2, anchor='N')
    pos1 = ax1.get_position()
    ax2.set_position([pos1.x0, 0.055, pos1.width, 0.025], which='both')
    bounds = np.arange(5,100,5)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ticks = np.arange(5,100,5)
    ticklabels = [str(int(x)) for x in ticks]
    cb1 = mpl.colorbar.ColorbarBase(
        ax2,
        cmap=cmap,
        norm=norm,
        boundaries=np.append(
            np.append([-1000.0],bounds), 
            [1000.0]
        ),
        extend='both',
        ticks=bounds,
        spacing='uniform',
        orientation='horizontal'
    )
    cb1.set_ticks(ticks)
    cb1.set_ticklabels(ticklabels)
    ax2.tick_params(axis='x', direction='in', width=1, length=6,
                    labelsize='small', bottom='on', top='off', pad=2)
    cb1.outline.set_linewidth(1.0)
    cb1.set_label('NMEP Reflectivity in range {args.lower_threshold}-{args.upper_threshold}dBZ',
                  fontweight='bold',fontsize='small',
                  horizontalalignment='right',x=1.0)

    plt.sca(ax1)

    CFT2 = ax1.contourf(
        datadict['meta']['xlons'], 
        datadict['meta']['xlats'], 
        datadict['data'], 
        bounds, 
        transform=datadict['meta']['projection'], 
        cmap=cmap, extend='both')

    valid_date = datadict['meta']['valid_date']
    plt.savefig(f"{figurename}_{datadict['meta']['valid_date']:%Y%m%dT%H%M}.png",
                dpi=150,
                facecolor='w',
                edgecolor='w',
                orientation='landscape',
                bbox_inches=None,
                pad_inches=0.5)
    plt.close('all')

def parse_args():
    parser = argparse.ArgumentParser(
        description = "plot the NMEP forecast ensemble for reflectivity"
    )
    parser.add_argument(
        "--wrf-grib-files",
        required=True,
        nargs="+",
        help="Absolute paths to wrf ensemble member gribs"
    )
    parser.add_argument(
        "--output-directory",
        default=os.getcwd(),
        help="Path to directory to output plots to"
    )
    parser.add_argument(
        '--name',
        help = 'A descriptive name for the run',
        default='UNKOWN'
    )
    parser.add_argument(
        '--lower-threshold',
        help='The lower dBZ threshold. Default is 30dBZ',
        default=30,
        type=int
    )
    parser.add_argument(
        '--upper-threshold',
        help='The upper dBZ threshold. If unspecified, defaults to 1000 which'
             ' effectively produces prob of exceedance (lower threshold) field.',
        default=1000,
        type=int
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()

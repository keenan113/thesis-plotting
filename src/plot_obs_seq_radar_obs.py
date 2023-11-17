import argparse
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from metpy.plots import ctables

import cartopy.feature as cfeature
import cartopy.crs as crs
import cartopy
from cartopy.feature import NaturalEarthFeature

import sys,os,shutil
from glob import glob
import warnings
warnings.filterwarnings("ignore")

import xarray
import pandas
from cartopy_util import get_cartopy_metadata

CopyIndices = {'observations':0,'prior_ens_mean':1,'posterior_ens_mean':2}

PlotStuff = {'refl':{   'VarName':'RADAR_REFLECTIVITY',
                        'cmap':ctables.registry.get_colortable('NWSReflectivity'),
                        'colorbar_title':'Reflectivity    $dB$',
                        'colorbar_interval':numpy.arange(5,80,5)
                    },
            'velo':{     'VarName':'DOPPLER_RADIAL_VELOCITY',
                        'cmap':plt.cm.bwr,
                        'colorbar_title':'Radial Velocity    $ms^{-1}$',
                        'colorbar_interval':numpy.arange(-60,64,4)
                    }
            }
BottomHeights = [250.,750.,1750.,2750.]
TopHeights = [750.,1250.,2250.,3250.]


def main():
    args = parse_args()

    variable_name = PlotStuff[args.field_name]['VarName']
    copy_index = CopyIndices[args.field_type]

    cartopy_metadata = get_cartopy_metadata(args.sample_file_path)
    """
    try:
        
        wrf_ds = xarray.open_dataset(args.sample_file_path)
        cart_lats = wrf_ds.XLAT.values[0,:,:]
        cart_lons = wrf_ds.XLONG.values[0,:,:]
        cart_proj = get_cartopy_crs(wrf_ds)
        xlim, ylim = calc_cart_extents(cart_proj,cart_lats,cart_lons)
        cartopy_metadata = cart_lats,cart_lons,cart_proj,xlim,ylim
    finally:
        print('Closing sample WRF file now...')
        wrf_ds.close()
    """
    for filename in args.obs_files:
        try:
            ds = read_obs_epoch_file(filename)
            variable_ds = subset_ds_to_variable(ds,variable_name)
        finally:
            ds.close()
        time = get_rounded_average_datetime(variable_ds.time)
        print("Working on {0:%Y%m%d_%H%M}".format(time))
        
        for height_min,height_max in zip(BottomHeights,TopHeights):
            print(f'\tFor {variable_name} {args.field_type} between {height_min} and {height_max}')
            height_ds = subset_ds_by_height_range(variable_ds,height_min,height_max)
            
            lon_da,lat_da,height_da = get_coordinates_from_ds(height_ds)
 
            plot_metadata = int(height_min),int(height_max),args.field_name,args.field_type,time
            
            plot_data = height_ds.isel(copy=copy_index).observations            

            if (~plot_data.notnull()).any():
                print('....skipping because we found nans')
                continue
            print(args.outdir)
            PlotReflObs(
                args.outdir,
                plot_data.values,
                lon_da.values,
                lat_da.values,
                plot_metadata,
                cartopy_metadata
            )   
def subset_ds_by_height_range(ds,height_min,height_max):
    height_ds = ds.where(
        (ds.location.isel(locdim=2) >= height_min) & (ds.location.isel(locdim=2) <= height_max),
        drop = True
    )
    return height_ds

def calc_cart_extents(dest_crs,lat,lon):
    # Need to modify the extents for the new projection
    bottom_left_lon = lon[0,0]
    bottom_left_lat = lat[0,0]
    top_right_lon = lon[-1,-1]
    top_right_lat = lat[-1,-1]
    src_crs = crs.PlateCarree()
    xs, ys, _ = dest_crs.transform_points(
        src_crs,
        numpy.array([bottom_left_lon, top_right_lon]),
        numpy.array([bottom_left_lat, top_right_lat])
    ).T
    _xlimits = xs.tolist()
    _ylimits = ys.tolist()

    return (_xlimits, _ylimits)

def get_cartopy_crs(wrfds):
    WRF_EARTH_RADIUS = 6370000.

    globe = crs.Globe(
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

def read_obs_epoch_file(filename):
    ds = xarray.open_dataset(filename)
    return ds

def subset_ds_to_variable(ds,variable_name):
    ob_type_indices = ds.where(
        ds.ObsTypesMetaData.astype('str').str.strip() == variable_name,
        drop=True
    ).ObsTypesMetaData.coords['ObsTypes'].values
    assert len(ob_type_indices)==1
    ob_type_index = ob_type_indices[0]
    ds = ds.where(ds.obs_type == ob_type_index,drop=True)
    return ds

def get_coordinates_from_ds(ds):
    #separate lon,lat, and height arrays so we can later subset data
    # arrays by height value

    location = ds.location

    lon = location.sel(locdim=0)
    lat = location.sel(locdim=1)
    height = location.sel(locdim=2)

    return lon,lat,height

def get_rounded_average_datetime(time_da):
    return pandas.to_datetime(time_da.mean().dt.round('1min').values)

def PlotReflObs(Outdir,ObsLevel,LonLevel,LatLevel,PlotSettings,cartopy_metadata):
    #lats, lons, projection, xlim, ylim = cartstuff
    HMin,HMax,VarKeyName,CopyKeyName,ValTime = PlotSettings
    cbarIntvl = PlotStuff[VarKeyName]['colorbar_interval']
    cbarTitle = PlotStuff[VarKeyName]['colorbar_title']
    VarName = PlotStuff[VarKeyName]['VarName']


    fig = plt.figure(figsize=(9, 10))
    ax1 = fig.add_subplot(2, 1, 1, projection=cartopy_metadata['projection'], anchor='S')
    ax1.set_position([0.05, 0.14, 0.9, 0.80], which='both')
    ax1.outline_patch.set_linewidth(1.0)
    ax1.outline_patch.set_zorder(6)
    ax1.set_title('Valid: {0:%m/%d/%Y @ %H:%M} UTC'.format(ValTime),fontsize='large', loc='right')
    ax1.set_title('{1} {0} b/w {2}m and {3}m'.format(CopyKeyName,VarName.lower(),HMin,HMax), fontsize='medium', loc='left')

    #Import the shapefiles for plotting the state lines, lakes, and coastlines
    states = NaturalEarthFeature(category='cultural', scale='10m',
                        facecolor='none', name='admin_1_states_provinces_shp')
    lakes = NaturalEarthFeature(category='physical', scale='10m',
                        facecolor='none', name='lakes', edgecolor='#ffffff')

    #Add the coastlines and states to the main axes
    ax1.add_feature(states, edgecolor='#404040', linewidth=2.0,
                        facecolor='none', zorder=4)
    ax1.add_feature(lakes, edgecolor='#404040', linewidth=2.0, zorder=4)

    #Set the x and y limits to equal the domain of the WRF run that we did
    ax1.set_xlim(cartopy_metadata['xlim'])
    ax1.set_ylim(cartopy_metadata['ylim'])

    cmap = PlotStuff[VarKeyName]['cmap']
    cmap.set_under('white')

    ax2 = fig.add_subplot(2, 1, 2 ,anchor='N')
    pos1 = ax1.get_position()
    ax2.set_position([pos1.x0, 0.055, pos1.width, 0.025], which='both')
    bounds = cbarIntvl
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ticks = cbarIntvl
    ticklabels = [str(int(x)) for x in ticks]
    cb1 = mpl.colorbar.ColorbarBase(ax2,
                                    cmap=cmap,
                                    norm=norm,
                                    boundaries=numpy.append(numpy.append([-100.0],
                                        bounds), [100.0]),
                                    extend='both',
                                    ticks=bounds,
                                    spacing='uniform',
                                    orientation='horizontal')
    cb1.set_ticks(ticks)
    cb1.set_ticklabels(ticklabels)
    ax2.tick_params(axis='x', direction='in', width=1, length=6,
                    labelsize='small', bottom='on', top='off', pad=2)
    cb1.outline.set_linewidth(1.0)
    cb1.set_label(cbarTitle,
                    fontweight='bold',
                    fontsize='small',
                    horizontalalignment='right',
                    x=1.0)
    #Call the main axes to the forefront
    plt.sca(ax1)

    #Plot the reflectivity colorfilled contours
    CFT2 = ax1.tricontourf(
        LonLevel,
        LatLevel,
        ObsLevel,
        bounds,
        cmap=cmap,
        transform=crs.PlateCarree(),
        extend='both'
    )

    figure_outpath = os.path.join(Outdir,'{4}_{3}_{0}m_to_{1}m_{2:%Y%m%d_%H%M}.png'.format(HMin,HMax,ValTime,CopyKeyName,VarName))
    print(figure_outpath)
    #Save the figure, the figure name is fed in from the function and the date labeling is handled here
    plt.savefig(figure_outpath)

    #Close all the plots that have been opened now that we have saved the plot
    plt.close('all')

def parse_args():
    parser = argparse.ArgumentParser(
        description = 'Plots the radar observaitons in observation space. Can plot observations, as well as ensemble estimates in observation space.'
    )
    parser.add_argument(
        '--obs-files',
        required=True,
        nargs = "+",
        help = 'List of obs_epoch_xxx.nc files. Typically produced with find'
    )
    parser.add_argument(
        '--outdir',
        default = os.path.join(os.getcwd(),"obsseq_radarplots"),
        help = "path to where to write plots"
    )
    parser.add_argument(
        "--field-name",
        required = True,
        choices = ["refl","velo"],
        help = "What type of radar observation to plot?"
    )
    parser.add_argument(
        "--field-type",
        default = "observations",
        choices = ["observations","prior_ens_mean","posterior_ens_mean"],
        help = "The type of plot to generate -- observations, model values, etc."
    )
    parser.add_argument(
        "--sample-file-path",
        required = True,
        help = "Path to any wrfout file that describes the domain to plot on"
    )
    return parser.parse_args()

if __name__ == '__main__':
    main()

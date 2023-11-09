import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from matplotlib.cm import get_cmap
import numpy as np

def plot_sea_level_pressure(figurename, title, datadict):

    fmtstr = '%3.0f'

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(2, 1, 1, projection=datadict['meta']['projection'], anchor='S')
    ax1.set_position([0.05, 0.09, 0.9, 0.9], which='both')
    #ax1.outline_patch.set_linewidth(1.0)
    #ax1.outline_patch.set_zorder(6)

    #Set plot extent using projection from src grid
    xlons = datadict['meta']['xlons']
    xlats = datadict['meta']['xlats']
    ax1.set_extent(
        [xlons[0,0],xlons[-1,-1],xlats[0,0],xlats[-1,-1]], 
        datadict['meta']['projection']
    )

    #Add state and lake borders to map
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

    #Set plot title(s) - with timing information on the right and actual
    # plot title on the left
    ax1.set_title(f"Valid: {datadict['meta']['valid_date']:%m/%d/%Y @ %H:%M} UTC",
                    fontsize='large',loc='right')
    ax1.set_title(title, fontsize='medium', loc='left')

    CFT2 = ax1.contour(
        datadict['meta']['xlons'], 
        datadict['meta']['xlats'], 
        datadict['data'], 
        np.arange(600,2000,10), 
        transform=crs.PlateCarree(), 
        colors='black',
        linewidths=0.75)
    print(f"{figurename}_{datadict['meta']['valid_date']:%Y%m%dT%H%M}.png")
    plt.savefig(f"{figurename}_{datadict['meta']['valid_date']:%Y%m%dT%H%M}.png",
                dpi=150,
                facecolor='w',
                edgecolor='w',
                orientation='landscape',
                bbox_inches=None,
                pad_inches=0.5)
    plt.close('all')
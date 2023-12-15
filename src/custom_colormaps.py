from metpy.plots import ctables
from matplotlib.cm import get_cmap
import matplotlib.colors as colors
import numpy

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=10):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

def get_reflectivity_colormap():
    cmap = ctables.registry.get_colortable('NWSReflectivity')
    cmap.set_under('white')
    return cmap

def get_fgen_colormap():
    cmap_tot = get_cmap('gnuplot2_r')
    cmap = truncate_colormap(cmap_tot, 0., 0.7)
    return cmap

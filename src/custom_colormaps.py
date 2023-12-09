from metpy.plots import ctables

def get_reflectivity_colormap():
    cmap = ctables.registry.get_colortable('NWSReflectivity')
    cmap.set_under('white')
    return cmap

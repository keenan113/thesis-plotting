import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from plot_cartopy import PlotCartopy

class PlotMap(PlotCartopy):
    
    def __init__(self):
        self.proj = None

    def plot_map(self,data,latitude,longitude,ax=None,save_path=None,prop={},map_prop={}):
        #Set default properties
        default_prop={'dots':True,'fillcolor':'category','cmap':None,'levels':None,'linecolor':'k','linewidth':1.0,'ms':7.5,'plot_names':False}
        default_map_prop={'res':'110m','land_color':'#FBF5EA','ocean_color':'#EDFBFF','linewidth':0.5,'linecolor':'k','figsize':(14,9),'dpi':200,'plot_gridlines':True}
        
        #Initialize plot
        prop = self.add_prop(prop,default_prop)
        map_prop = self.add_prop(map_prop,default_map_prop)
        self.plot_init(ax,map_prop)

        #try:
        bounds = self.set_projection(latitude,longitude)
        #except:
        #    print('No lat/lon information found. Defaulting to global plot')
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt


class PlotCartopy:

    def check_res(self,res):
        if res not in ['110m','50m','10m']:
            print('Expected a resolution string like 10m, 50m or 110m. Saw none so defaulting to 110m')
            return '110m'
        else:
            return res

    def create_cartopy(self,proj='PlateCarree',mapobj=None,**kwargs):
    
        """
        Initialize a cartopy instance passed projection.
        
        Parameters:
        -----------
        projection
            String representing the cartopy map projection.
        ax
            Axis on which to draw on. Default is None.
        mapobj
            Existing cartopy projection. If passed, will be used instead of generating a new one.
        **kwargs
            Additional arguments that are passed to those associated with projection.
        """
        
        #Initialize an instance of cartopy if not passed
        if mapobj is None:
            self.proj = getattr(ccrs, proj)(**kwargs)
        else:
            self.proj = mapobj

    def create_geography(self,prop):
        
        """
        Set up the map geography and colors.
        
        Parameters:
        -----------
        prop : dict
            dict entry containing information about the map geography and colors
        """
        
        #get resolution corresponding to string in prop
        res = self.check_res(prop['res'])
        
        #Add "hidden" state alpha prop
        if 'state_alpha' not in prop.keys(): prop['state_alpha'] = 1.0
        
        #fill oceans if specified
        self.ax.set_facecolor(prop['ocean_color'])
        ocean_mask = self.ax.add_feature(cfeature.OCEAN.with_scale(res),facecolor=prop['ocean_color'],edgecolor='face')
        lake_mask = self.ax.add_feature(cfeature.LAKES.with_scale(res),facecolor=prop['ocean_color'],edgecolor='face')
        continent_mask = self.ax.add_feature(cfeature.LAND.with_scale(res),facecolor=prop['land_color'],edgecolor='face')
        
        #draw geography
        states = self.ax.add_feature(cfeature.STATES.with_scale(res),linewidths=prop['linewidth'],linestyle='solid',edgecolor=prop['linecolor'],alpha=prop['state_alpha'])
        countries = self.ax.add_feature(cfeature.BORDERS.with_scale(res),linewidths=prop['linewidth'],linestyle='solid',edgecolor=prop['linecolor'])
        coastlines = self.ax.add_feature(cfeature.COASTLINE.with_scale(res),linewidths=prop['linewidth'],linestyle='solid',edgecolor=prop['linecolor'])

    def plot_init(self,ax,map_prop,plot_geography=True):
        
        """
        Initializes the plot by creating a cartopy and axes instance, if one hasn't been created yet, and adds geography.
        
        Parameters:
        -----------
        ax : axes
            Instance of axes
        map_prop : dict
            Dictionary of map properties
        """

        #create cartopy projection, if none existing
        if self.proj is None:
            self.create_cartopy(proj='PlateCarree',central_longitude=0.0)
        
        #create figure
        if ax is None:
            self.fig = plt.figure(figsize=map_prop['figsize'],dpi=map_prop['dpi'])
            self.ax = plt.axes(projection=self.proj)
        else:
            fig_numbers = [x.num for x in mlib._pylab_helpers.Gcf.get_all_fig_managers()]
            if len(fig_numbers) > 0:
                self.fig = plt.figure(fig_numbers[-1])
            else:
                self.fig = plt.figure(figsize=map_prop['figsize'],dpi=map_prop['dpi'])
            self.ax = ax
        
        #Attach geography to plot, lat/lon lines, etc.
        if plot_geography:
            self.create_geography(map_prop)

    def add_prop(self,input_prop,default_prop):
        
        r"""
        Overrides default property dictionary elements with those passed as input arguments.
        
        Parameters:
        -----------
        input_prop : dict
            Dictionary to use for overriding default entries.
        default_prop : dict
            Dictionary containing default entries.
        
        Returns:
        --------
        dict
            Default dictionary overriden by entries in input_prop.
        """
        
        #add kwargs to prop and map_prop
        for key in input_prop.keys(): default_prop[key] = input_prop[key]
            
        #Return prop
        return default_prop

    def set_projection(self,latitude,longitude):
        bound_w = longitude[0,0]
        bound_e = longitude[-1,-1]
        bound_s = latitude[0,0]
        bound_n = latitude[-1,-1]
        #Set map extent
        self.ax.set_extent([bound_w,bound_e,bound_s,bound_n], crs=ccrs.PlateCarree())
            
        return bound_w, bound_e, bound_s, bound_n



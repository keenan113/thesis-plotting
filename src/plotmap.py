import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cartopy.feature import NaturalEarthFeature
import cartopy.crs as crs

import custom_colormaps

class PlotMap:
    def __init__(self, outdir, filename_prefix, cartopy_metadata, latitude, longitude):
   
        self.filename = None 
        self.outdir = outdir
        self.filename_prefix = filename_prefix
        self.plot_title_prefix = ' '.join(filename_prefix.split('_'))

        self.projection = cartopy_metadata.projection
        self.geobounds = cartopy_metadata.geobounds
        self.xlim = cartopy_metadata.xlim
        self.ylim = cartopy_metadata.ylim
        
        self.latitude = latitude
        self.longitude = longitude

        self.states = NaturalEarthFeature(category='cultural', scale='10m', facecolor='none',
                                          name='admin_1_states_provinces_shp')
        self.lakes = NaturalEarthFeature(category='physical', scale='10m', facecolor='none', name='lakes',
                                         edgecolor='#ffffff')
        self.contour_line_width = 2

        
    def init_plot(self):
        self.fig = plt.figure(figsize=(10, 10))
        self.ax1 = self.fig.add_subplot(2, 1, 1, projection=self.projection, anchor='S')
        self.ax2 = self.fig.add_subplot(2, 1, 2, anchor='N')
       
        self.ax1.set_position([0.05,0.14,0.9,0.80], which='both') 
        self.ax1.outline_patch.set_linewidth(1.0)
        self.ax1.outline_patch.set_zorder(6)
        
        self._add_map_features()

        self.ax1.set_extent(
            self._define_plot_extents(),
            self.projection
        )
        self._set_xylims()
        
    def _set_xylims(self):
        self.ax1.set_xlim(self.xlim)
        self.ax1.set_ylim(self.ylim)
    
    def _define_plot_extents(self):
        plot_extents = [
            self.geobounds.bottom_left.longitude,
            self.geobounds.top_right.longitude,
            self.geobounds.bottom_left.latitude,
            self.geobounds.top_right.latitude
        ]
        return plot_extents

    def _add_map_features(self):
        self.ax1.add_feature(self.states, edgecolor='#404040', linewidth=1.0, facecolor='none', zorder=4)
        self.ax1.add_feature(self.lakes, edgecolor='#404040', linewidth=1.0, zorder=4)
    
    def _define_colormap(self,cmap_name):
        if cmap_name == 'reflectivity':
            self.cmap = custom_colormaps.get_reflectivity_colormap()
        elif cmap_name == "fgen":
            self.cmap = custom_colormaps.get_fgen_colormap()
        else:
            self.cmap = plt.get_cmap('viridis')
    
    def _configure_colormap_legend(self,bounds=np.arange(1,1000,100),ticks=np.arange(1,1000,100),extend='both',legend_name="contour field"):

        pos1 = self.ax1.get_position()

        self.ax2.set_position([pos1.x0, 0.055, pos1.width, 0.025], which='both')

        norm = mpl.colors.BoundaryNorm(bounds, self.cmap.N)
        ticklabels = [str(int(x)) for x in ticks]
        cb1 = mpl.colorbar.ColorbarBase(
            self.ax2,
            cmap=self.cmap,
            norm=norm,
            boundaries=np.append(np.append([-100.0], bounds), [100.0]),
            extend=extend,
            ticks=bounds,
            spacing='uniform',
            orientation='horizontal'
        )
        cb1.set_ticks(ticks)
        cb1.set_ticklabels(ticklabels)
        self.ax2.tick_params(
            axis='x',
            direction='in',
            width=1,
            length=6,
            labelsize='small',
            bottom='on',
            top='off',
            pad=2
        )
        cb1.outline.set_linewidth(1.0)
        cb1.set_label(
            legend_name,
            fontweight='bold',
            fontsize='small',
            horizontalalignment='right',
            x=1.0
        )

    def _set_plot_custom_title(self,figure_title=None):
        figure_titles = []
        if self.plot_title_prefix:
            figure_titles.append(self.plot_title_prefix)

        if figure_title is not None:
            figure_titles.append(figure_title)
        
        self.ax1.set_title('\n'.join(figure_titles), fontsize='medium', loc='left')


    def add_filled_contour(self,contour_data,contour_interval,cmap):
        
        filled_contour = self.ax1.contourf(
            self.longitude,
            self.latitude,
            contour_data,
            contour_interval,
            transform=crs.PlateCarree(),
            cmap=cmap,
            extend='both'
        )
        return filled_contour
     
    def add_contour(self,contour_data,contour_interval,contour_color):
        contour = self.ax1.contour(
            self.longitude,
            self.latitude,
            contour_data,
            contour_interval,
            linewidths = self.contour_line_width,
            colors = contour_color,
            transform = crs.PlateCarree()
        )
        return contour

    def _set_plot_time_title(self,valid_time,issuance_time=None):
        title_list = [f"Valid: {valid_time:%m/%m/%Y %H:%M} UTC"]
        if issuance_time:
            title_list.append(f"Issued: {issuance_time:%m/%m/%Y %H:%M} UTC")
        title_str = '\n'.join(title_list)
        self.ax1.set_title(title_str, fontsize='large', loc='right')

    def set_contour_labels(self,contour_field):
        self.ax1.clabel(contour_field, inline=1, fontsize='small', fmt='%3.0f')

    def save_plot(self):
        if self.filename == None:
            self.filename = "./test.png"
        plt.savefig(self.filename, dpi=150, facecolor='w', edgecolor='w',
                    orientation='landscape', bbox_inches=None, pad_inches=0.5)
        plt.close('all')
    
    def plot_MSLP_GH500mb_CREFL(self,valid_time,MSLP,GH500mb,CREFL,issuance_time=None):
         
        self.init_plot()
        self._set_plot_time_title(valid_time,issuance_time=issuance_time)
        self._set_plot_custom_title(figure_title=f'CREFL [color-filled] | MSLP [black] | 500hPa GH [red]') 
        
        self._define_colormap('reflectivity') 
        self._configure_colormap_legend(
            legend_name='Reflectivity    $dB$',
            bounds=np.arange(5, 80, 5),
            ticks=np.arange(5, 80, 5)
        )
            

        MSLP_contours = self.add_contour(MSLP, np.arange(980, 1040, 2),'black')
        GH500mb_contours = self.add_contour(GH500mb, np.arange(4800, 6000, 100), 'red')
        CREFL_contourfilled = self.add_filled_contour(CREFL, np.arange(5, 80, 5), self.cmap)

        self.set_contour_labels(MSLP_contours)
        self.set_contour_labels(GH500mb_contours)
        
        self._set_xylims()
        self._set_plot_filename("CREFL_MSLP_500mbGPT",valid_time,issuance_time=issuance_time)
        self.save_plot()

    def plot_isobaric_fgen_heights(self,valid_time,FGEN,GH,issuance_time=None):
         
        self.init_plot()
        self._set_plot_time_title(valid_time,issuance_time=issuance_time)
        self._set_plot_custom_title(figure_title=f'FGEN [color-filled] | {pressure_level}hPa GH [black]') 
        
        self._define_colormap('fgen') 
        self._configure_colormap_legend(
            legend_name='Frontogenesis (K/100km/3hr)',
            bounds=np.arange(1, 19, 2),
            ticks=np.arange(1, 19, 2)
        )
            
        
        GH_contours = self.add_contour(GH, np.arange(1000, 10000, 100),'black')
        FGEN_contourfilled = self.add_filled_contour(FGEN, np.arange(1, 19, 2), self.cmap)

        self.set_contour_labels(GH_contours)
        
        self._set_xylims()
        self._set_plot_filename(f"FGEN_GH_{pressure_level}",valid_time,issuance_time=issuance_time)
        self.save_plot()


    def _set_plot_filename(self,figure_name,valid_time,issuance_time=None):
        valid_time_chunk = f"_valid_{valid_time:%Y%m%dT%H%M}"
        if issuance_time is None:
            issuance_time_chunk = ""
        else:
            issuance_time_chunk = f"_issued_{issuance_time:%Y%m%dT%H%M}"
        self.filename  = os.path.join(
            self.outdir,
            f"{self.filename_prefix}_{figure_name}{issuance_time_chunk}{valid_time_chunk}.png"
        )




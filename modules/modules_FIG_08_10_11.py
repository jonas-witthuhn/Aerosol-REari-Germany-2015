import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def make_map(fig,ax):
    def add_features(ax):
        """ Add Map features for a nice looking map. -> Country Boundaries, Lakes, Rivers"""
        countries = cfeature.NaturalEarthFeature(category='cultural',
                                                        name='admin_0_countries',
                                                       scale='50m',
                                                       facecolor='none')
        lakes = cfeature.NaturalEarthFeature(category='physical',
                                             name='lakes',
                                            scale='10m')
        lakes_euro = cfeature.NaturalEarthFeature(category='physical',
                                             name='lakes_europe',
                                            scale='10m')

        rivers = cfeature.NaturalEarthFeature(category='physical',
                                             name='rivers_lake_centerlines',
                                            scale='10m',
                                             facecolor='none')

        rivers_euro = cfeature.NaturalEarthFeature(category='physical',
                                             name='rivers_europe',
                                            scale='10m',
                                             facecolor='none')

        ax.add_feature(countries,edgecolor='k')
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
        ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['water'])

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        gl.top_labels = False
        gl.left_labels = True
        gl.right_labels = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.xformatter = LONGITUDE_FORMATTER
        return ax

    ax.gridlines()
    extent = [5.8,15.2,47.2,55]
    ax.set_extent(extent,ccrs.PlateCarree())
    ax = add_features(ax)
    return fig,ax

def make_plot(fig,sp,x,y,data,data2=None,title='',
              clabel='',cmap='autumn_r',cbaxextend=[],
              cbpanchor=False,cbspacing='proportional',
              cb=True,cbarres=2,cbticks=None,cbaspect=50,
              data2res=2,data2ticks=None):
    # make subplot
    ax = fig.add_subplot(*sp,projection=ccrs.PlateCarree())
    ax.set_title(title,fontsize=14)
    
    # plot map
    fig,ax = make_map(fig,ax)
   
    if type(cbticks)==type(None):
        # make contourplot
        ticks = np.arange(np.round(np.nanpercentile(data,1),cbarres)-10**(-cbarres),
                          np.round(np.nanpercentile(data,99),cbarres)+10**(-cbarres),10**(-cbarres))
        while len(ticks)>11:
            ticks = ticks[::2]
    else:
        ticks=cbticks
    p1=ax.contourf(x,y,data,transform=ccrs.PlateCarree(),cmap=cmap,
                   levels = ticks)
    
    if type(data2)!=type(None):
        if type(data2ticks)==type(None):
            # make contourplot
            ticks2 = np.arange(np.round(np.nanpercentile(data2,1),data2res),
                               np.round(np.nanpercentile(data2,99),data2res),10**(-data2res))
            while len(ticks2)>11:
                ticks2 = ticks2[::2]
        else:
            ticks2=data2ticks
        
        CS=ax.contour(x,y,data2,transform=ccrs.PlateCarree(),levels=ticks2,colors=['k']*len(ticks2))
        ax.clabel(CS, inline=1, fontsize=14)
    
    if cb:
        cbaxes = ax
        if len(cbaxextend)>0:
            cbaxes = [ax]+cbaxextend
        if type(cbpanchor)==type(bool):
            cb = plt.colorbar(p1,ax=cbaxes,ticks=ticks,
                              aspect=cbaspect,
                              spacing=cbspacing,pad=0.03,shrink=0.9,
                              label = ' ')
        else:
            cb = plt.colorbar(p1,ax=cbaxes,ticks=ticks,
                              panchor=cbpanchor,
                              aspect=cbaspect,
                              spacing=cbspacing,pad=0.03,shrink=0.9,
                              label = ' ')
        cb.ax.set_title(clabel,loc='left')
    fig.canvas.draw()
    return fig, ax
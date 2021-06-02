import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from matplotlib.patches import Rectangle

def get_metrics(O,P):
    #ensure no nans in data
    idx = ~np.isnan(O)*~np.isinf(O)
    idx*= ~np.isnan(P)*~np.isinf(P)
    O = O[idx]
    P = P[idx]
    # precalculation
    N = len(O)
    N1 = 1./float(N)
    N2 = 1./(float(N)+1.)
    DELTA = P-O
    SUM = P+O
    # MBE
    MBE = N1 * np.sum(DELTA)
    # MABE
    MABE = N1 * np.sum(np.abs(DELTA))
    # RMSE
    RMSE = np.sqrt(N2*np.sum(DELTA**2))
    # fractional bias FB or normalzed mean bias MNMB
    FB = 2.*N1 * np.sum(DELTA / SUM)
    # fractional gross error FGE
    FGE = 2.*N1 * np.sum(np.abs(DELTA)/np.abs(SUM))
    return MBE,MABE,RMSE,FB,FGE

def make_plotqq(fig,
              X,Y,W=None,
              ax=None,
              ax2=None,
              lalpha=1,
              title='',title2='',
              Xlabel='',Ylabel='',Wlabel='',Value='',
              xlim=(0,1),xlim2=(0,1),ylim=(0,1),
              qqplot=True,percentiles=None,
              histplot=True,metrics_loc=0):

    idx = ~np.isnan(X)*~np.isinf(X)
    idx*= ~np.isnan(Y)*~np.isinf(Y)
    X = X[idx]
    Y = Y[idx]
    
    if type(ax)!=type(None) and not qqplot:
        ax.set_title(title,fontsize=14)
        if type(percentiles)==None:
            percentiles = np.arange(1,101)
        else:
            percentiles = np.sort(np.array(percentiles))
        x = np.percentile(stats.norm.rvs(loc=0, scale=1, size=1000000),percentiles)
        Xstd = np.nanstd(X)
        Ystd = np.nanstd(Y)
        # sort by X
        ind = np.argsort(X)
        pY = Y[ind]
        pX = X[ind]
        # percentile of sorted X point
        Xpercentiles = 100* np.arange(len(X))/(len(X)-1)
        Xnormp = np.percentile(stats.norm.rvs(loc=0, scale=1, size=1000000),Xpercentiles)

        # look for X values at special Xnormp
        ispecial = np.searchsorted(Xnormp,[-2,-1,0,1,2])

        # mean delta at percentile step
        DELTA = pY-pX
        pDELTA = np.zeros(len(percentiles)-1)
        for i,p in enumerate(percentiles[1:]):
            idx = Xpercentiles<=p
            idx*= Xpercentiles>percentiles[i]
            pDELTA[i] = np.nanmean(DELTA[idx])

        # calculate data percentiles
        xp = np.nanpercentile(X,percentiles)
        xp50 = np.nanpercentile(X,50)
        yp = np.nanpercentile(Y,percentiles)
        yp50 = np.nanpercentile(Y,50)

        # colormesh 
        gridx=np.linspace(xlim[0],xlim[1],50)
        gridy=np.linspace(ylim[0],ylim[1],50)
        grid,_,_=np.histogram2d(Xnormp,DELTA, bins=[gridx, gridy]) # number of points in bins
        grid[grid==0]=np.nan
        p1=ax.pcolormesh(gridx,gridy,grid.T,cmap='gnuplot',alpha=0.5,zorder=1)
        cbar = fig.colorbar(p1,ax=ax)
        cbar.ax.set_title(Wlabel,fontsize=14)
        # calculate bin mean
        wgrid,_,_=np.histogram2d(Xnormp,DELTA, bins=[gridx, gridy],weights=DELTA)
        Wgrid = np.nansum(wgrid,axis=1)
        Ngrid = np.nansum(grid,axis=1)
        ax.plot(gridx[:-1]+0.5*(gridx[1:]-gridx[:-1]),Wgrid/Ngrid,color='k',linewidth=2,zorder=2)

        # 0 line
        ax.axhline(0,color='k',linestyle=':',linewidth=2,zorder=3)

        # X values at full quantiles:
        xspan = xlim[1]-xlim[0]
        yspan = ylim[1]-ylim[0]
        ax.add_patch(Rectangle((xlim[0],ylim[0]),xspan,0.19*yspan,
                                fill=True,facecolor='w',edgecolor='k',alpha=0.7,zorder=3))
        ax.annotate(f"{Xlabel} {Value}:",(xlim[0]+0.01*xspan,ylim[0]+0.09*yspan),ha='left',va='bottom',fontsize=12)#,fontweight='bold'
        for i in ispecial:
            ax.annotate(f"{X[ind][i]:.2f}",(Xnormp[i],ylim[0]+0.005*yspan),ha='center',va='bottom',fontsize=12)#,fontweight='bold'


        ax.set_xlabel(r'$\sigma$'+f' from {Xlabel} median [-]',fontsize=14)
        ax.set_ylabel(r'$\Delta\,$'+f'{Value}',fontsize=14)
        ax.set_ylim(ylim)
        ax.grid(True)
        
    if qqplot:
        ax.set_title(title,fontsize=14)
        if type(percentiles)==None:
            percentiles = np.arange(1,101)
        else:
            percentiles = np.sort(np.array(percentiles))
        x = np.percentile(stats.norm.rvs(loc=0, scale=1, size=1000000),percentiles)
        Xstd = np.nanstd(X)
        Ystd = np.nanstd(Y)
        # sort by X
        ind = np.argsort(X)
        pY = Y[ind]
        pX = X[ind]
        # percentile of sorted X point
        Xpercentiles = 100* np.arange(len(X))/(len(X)-1)
        Xnormp = np.percentile(stats.norm.rvs(loc=0, scale=1, size=1000000),Xpercentiles)
        
        # look for X values at special Xnormp
        ispecial = np.searchsorted(Xnormp,[-2,-1,0,1,2])
        
        # mean delta at percentile step
        DELTA = pY-pX
        pDELTA = np.zeros(len(percentiles)-1)
        for i,p in enumerate(percentiles[1:]):
            idx = Xpercentiles<=p
            idx*= Xpercentiles>percentiles[i]
            pDELTA[i] = np.nanmean(DELTA[idx])
        
        # calculate data percentiles
        xp = np.nanpercentile(X,percentiles)
        xp50 = np.nanpercentile(X,50)
        yp = np.nanpercentile(Y,percentiles)
        yp50 = np.nanpercentile(Y,50)
        
        # colormesh 
        gridx=np.linspace(xlim[0],xlim[1],50)
        gridy=np.linspace(ylim[0],ylim[1],50)
        grid,_,_=np.histogram2d(Xnormp,DELTA, bins=[gridx, gridy]) # number of points in bins
        grid[grid==0]=np.nan
        p1=ax.pcolormesh(gridx,gridy,grid.T,cmap='gnuplot',alpha=0.5,zorder=1)
        
        # mean DELTA
        ax.plot(x[1:],pDELTA,color='k',linewidth=2,zorder=2)
        ax.axhline(0,color='k',linewidth=3,zorder=3)
        ax.set_ylabel(r'$\Delta\,$'+f'{Value}',fontsize=14)
        
        # QQ-plot
        ax1 = ax.twinx()
        xstd = (xp-xp50)/Xstd
        ystd = (yp-yp50)/Ystd

        ax1.plot(xstd[1:-1],ystd[1:-1],label='Quantile relation',linewidth=3,zorder=2)
        ax1.plot(xlim,xlim,'k',alpha=0.5,linewidth=2,linestyle='--',zorder=2)
        ax1.plot([-4,-4],[-4,-4],'k',label=r'$\Delta\,$'+f'{Value}',linewidth=2,linestyle='-',zorder=2)
        
        # X values at full quantiles:
        ax1.add_patch(Rectangle((-3,-3),6,0.72,fill=True,facecolor='w',edgecolor='k',alpha=0.5,zorder=3))
        ax1.annotate(f"AERONET {Value}:",(-2.9,-2.65),ha='left',va='bottom',fontsize=12,fontweight='bold')
        for i in ispecial:
            ax1.annotate(f"{X[ind][i]:.2f}",(Xnormp[i],-2.99),ha='center',va='bottom',fontsize=12,fontweight='bold')
        
        ax1.grid(True)
        ax.set_xlabel(r'AERONET quantiles [$\sigma$]',fontsize=14)
        ax1.set_ylabel(r'CAMS RA quantiles [$\sigma$]',fontsize=14)
        ax1.legend(fontsize = 14)
        ax1.set_xlim(xlim)
        ax1.set_ylim(xlim)
        ax.set_ylim(ylim)
        ax.grid(True)

    if histplot:
        ax2.set_title(title2,fontsize=14)
        ax2.hist((X,Y),20,range=xlim2,
                histtype='bar',
                log=False,
                label=(Xlabel,Ylabel))
        ax2.legend(loc=1,fontsize=14,framealpha=1,edgecolor='k')
        corr=np.corrcoef(X,Y)[0,1]
        N= len(Y)
        MBE,MABE,RMSE,FB,FGE = get_metrics(X,Y)
        met = str(f'Performance metrics:\n'+
                  f'  N    = {N}\n'+
                  f'  R    = {corr:.2f}\n'+
                  f'  MBE  = {MBE:.2f}\n'+
                  f'  RMSE = {RMSE:.2f}')

        if metrics_loc == 0:
            ax2.text(0.03,0.97,met,
                     fontsize=12,
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    bbox=dict(facecolor='w',alpha=lalpha))
        else:
            ax2.text(0.97,0.03,met,
                     fontsize=12,
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    multialignment='left',
                    transform=ax2.transAxes,
                    bbox=dict(facecolor='w',alpha=lalpha))          
        
        ax2.set_xlabel(Value,fontsize=14)
        ax2.set_ylabel('Count [-]',fontsize=14)
        ax2.grid(True)
    else:
        fig2 = None
        ax2 = None
    return fig      

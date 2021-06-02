import os
import xarray as xr
import numpy as np
import pandas as pd

import modules.load_data as ld



def sel(ds,st=None,seas=None,tag=None,sza=90.):
    """ selector function. select station, season, location-tagg and mask sza"""
    # select according sza
    ds = ds.where(ds.mu0>np.cos(np.deg2rad(sza)),drop=True)
    if len(ds.col)==0:
        return False
    if st != None:
        ds = ds.where(ds.name==st,drop=True)
    if len(ds.col)==0:
        return False
    if tag != None:
        ds = ds.sel(tag=tag)
        ds = ds.where(ds.mask,drop=True)
    if len(ds.col)==0:
        return False
    if seas != None:
        grp = ds.groupby('time.season').groups
        if not seas in grp.keys():
            return False
        ds = ds.isel(dict(col=grp[seas]))
    if len(ds.col)==0:
        return False
    return ds



def meanF(ds,typ='DNI',st=None,seas=None,tag=None,sza=90):
    """ calculate annual mean irradiance of input type"""
    ds = sel(ds,st,seas,tag,sza)
    if type(ds) == bool:
        return np.nan,np.nan,np.nan,np.nan
    else:
        if typ == 'DNI':
            N = np.count_nonzero(ds.dni_obs)
            Fmean = [N,
                     np.nanmean(ds.dni_obs),
                     np.nanmean(ds.dni_clr),
                     np.nanmean(ds.dni_aer)]
            return Fmean
        elif typ == 'GHI':
            N = np.count_nonzero(ds.glo_obs)
            Fmean = [N,
                     np.nanmean(ds.glo_obs),
                     np.nanmean(ds.glo_clr),
                     np.nanmean(ds.glo_aer)]
            return Fmean
        elif typ == 'DHI':
            N = np.count_nonzero(ds.glo_obs-(ds.dni_obs*ds.mu0))
            Fmean = [N,
                     np.nanmean(ds.glo_obs-(ds.dni_obs*ds.mu0)),
                     np.nanmean(ds.glo_clr-(ds.dni_clr*ds.mu0)),
                     np.nanmean(ds.glo_aer-(ds.dni_aer*ds.mu0))]
            return Fmean
        else:
            raise ValueError("Keyword 'typ' must be one of ['DNI','GHI','DHI']!")

def dsigma(ds,typ='DNI',st=None,seas=None,tag=None,sza=90):
    """ skill test """
    ds = sel(ds,st,seas,tag,sza)
    if type(ds) == bool:
        return np.nan,np.nan,np.nan
    else:
        if typ =='DNI':
            dclr_obs = ds.dni_clr-ds.dni_obs
            daer_obs = ds.dni_aer-ds.dni_obs
            dsigm_clr = float(np.std(dclr_obs))
            dsigm_aer = float(np.std(daer_obs))
            dsigm = dsigm_clr - dsigm_aer
            return dsigm,dsigm_clr,dsigm_aer
        elif typ =='GHI':
            dclr_obs = ds.glo_clr-ds.glo_obs
            daer_obs = ds.glo_aer-ds.glo_obs
            dsigm_clr = float(np.std(dclr_obs))
            dsigm_aer = float(np.std(daer_obs))
            dsigm = dsigm_clr - dsigm_aer
            return dsigm,dsigm_clr,dsigm_aer
        elif typ =='DHI':
            dhi_obs = ds.glo_obs-(ds.dni_obs*ds.mu0)
            dhi_clr = ds.glo_clr-(ds.dni_clr*ds.mu0)
            dhi_aer = ds.glo_aer-(ds.dni_aer*ds.mu0)
            dclr_obs = dhi_clr-dhi_obs
            daer_obs = dhi_aer-dhi_obs
            dsigm_clr = float(np.std(dclr_obs))
            dsigm_aer = float(np.std(daer_obs))
            dsigm = dsigm_clr - dsigm_aer
            return dsigm,dsigm_clr,dsigm_aer
        else:
            raise ValueError("Keyword 'typ' must be one of ['DNI','GHI','DHI']!")



def get_metrics(ds,typ,st=None,seas=None,tag=None,sza=90):
    """ calculate statistical metrics: correlation, standard deviation
        mean bias, RMSE, fractional bias and error. Calculate direct/diffuse 
        ratio"""
    ds = sel(ds,st,seas,tag,sza)
    if type(ds) == bool:
        return [np.nan]*6

    if typ == 'DNI':
        O = ds.dni_obs.values
        P = ds.dni_aer.values
    elif typ == 'GHI':
        O = ds.glo_obs.values
        P = ds.glo_aer.values
    elif typ == 'DHI':
        O = ds.glo_obs-(ds.dni_obs*ds.mu0)
        P = ds.glo_aer-(ds.dni_aer*ds.mu0)
        
        O=O.values
        P=P.values
    else:
        raise ValueError("Keyword 'typ' must be one of ['DNI','GHI','DHI']!") 

    idx = ~np.isnan(O)*~np.isinf(O)
    idx*= ~np.isnan(P)*~np.isinf(P)
    O = O[idx]
    P = P[idx]
    
    Dobs = np.nanmean(O)
    Daer = np.nanmean(P)
    
    # precalculation
    N = len(O)
    if N<3:
        return [np.nan]*6
    
    N1 = 1./float(N)
    DELTA = P-O
    SUM = P+O


    # Pearson correlation
    R = float(np.corrcoef(np.vstack((O,P)))[0,1])
    # MBE
    MBE = N1 * np.sum(DELTA)
    # RMSE
    RMSE = np.sqrt(N1*np.sum(DELTA**2))
    # fractional bias FB or normalzed mean bias MNMB
    FB = 2.*N1 * np.sum(DELTA / SUM)
    # fractional gross error FGE
    FGE = 2.*N1 * np.sum(np.abs(DELTA)/np.abs(SUM))

    return Dobs,Daer,N,R,MBE,RMSE,FB,FGE

def get_metrics_minus_clr(ds,typ,st=None,seas=None,tag=None,sza=90):
    """ calculate statistical metrics: correlation, standard deviation
    mean bias, RMSE, fractional bias and error. Calculate direct/diffuse
    ratio.
    Calculation are done on attenuated irradiance 
    (F with aerosol) - (F without aerosol)
    """
    ds = sel(ds,st,seas,tag,sza)
    if type(ds) == bool:
        return [np.nan]*6

    if typ == 'DNI':
        O = ds.dni_obs.values
        P = ds.dni_aer.values
        C = ds.dni_clr.values
    elif typ == 'GHI':
        O = ds.glo_obs.values
        P = ds.glo_aer.values
        C = ds.glo_clr.values
    elif typ == 'DHI':
        O = ds.glo_obs-(ds.dni_obs*ds.mu0)
        P = ds.glo_aer-(ds.dni_aer*ds.mu0)
        C = ds.glo_clr-(ds.dni_clr.values*ds.mu0)
        
        O=O.values
        P=P.values
        C=C.values
    else:
        raise ValueError("Keyword 'typ' must be one of ['DNI','GHI','DHI']!") 

    O = O-C
    P = P-C
        
    idx = ~np.isnan(O)*~np.isinf(O)
    idx*= ~np.isnan(P)*~np.isinf(P)
    O = O[idx]
    P = P[idx]
    
    # precalculation
    N = len(O)
    if N<3:
        return [np.nan]*6
    
    Dobs = np.nanmean(O)
    Daer = np.nanmean(P)
    
    N1 = 1./float(N)
    DELTA = P-O
    SUM = P+O

    # Pearson correlation
    R = float(np.corrcoef(np.vstack((O,P)))[0,1])
    # MBE
    MBE = N1 * np.sum(DELTA)
    # RMSE
    RMSE = np.sqrt(N1*np.sum(DELTA**2))
    # fractional bias FB or normalzed mean bias MNMB
    FB = 2.*N1 * np.sum(DELTA / SUM)
    # fractional gross error FGE
    FGE = 2.*N1 * np.sum(np.abs(DELTA)/np.abs(SUM))

    return Dobs,Daer,N,R,MBE,RMSE,FB,FGE
        
    


def metrics_table(typ,stations,dsy,sza=90,difdir=False):
    """ print Tables 04 05 06"""
    def _add_data(dsy,typ,st,seas,tag,sza,difdir=False):
        O,P,N,R,MBE,RMSE,FB,FGE = get_metrics(dsy,typ=typ,st=st,seas=seas,tag=tag,sza=sza)
        O0,P0,N0,R0,MBE0,RMSE0,FB0,FGE0 = get_metrics_minus_clr(dsy,typ=typ,st=st,seas=seas,tag=tag,sza=sza)
        dsig,sigclr,sigaer = dsigma(dsy,typ=typ,tag=tag,seas=seas,sza=sza)
        _,Fobs,Fclr,Faer = meanF(dsy,typ=typ,tag=tag,seas=seas,sza=sza)
        
        ## Line
        # --------|----------------------------------|---Faer vs Fobs--|--Faer vs obs-clr -| (sigma (clr-obs) - sigma(aer-obs))
        # Tag & N & meanFclr & mean Fobs & mean Faer &  R & RMSE & MBE & R & RMSE & MBE    & sigma
        line = ''
        line+= f' {N:d} & {Fclr:.0f} & {Fobs:.0f} & {Faer:.0f} &'
        if difdir:
            _,GLOobs,_,GLOaer = meanF(dsy,typ='GHI',tag=tag,seas=seas,sza=sza)
            _,DHIobs,_,DHIaer = meanF(dsy,typ='DHI',tag=tag,seas=seas,sza=sza)
            DDRobs = DHIobs / (GLOobs-DHIobs)
            DDRaer = DHIaer / (GLOaer-DHIaer)
            line+= f' {DDRobs:.3f} & {DDRaer:.3f} &'  
            
        line+= f' {RMSE:.0f} & {MBE:.0f} &'
        line+= f' {R:.3f} & {R0:.3f} &'
        line+= f' {dsig:.2f} \\\\'
        return line
    
    
    tags = xr.DataArray(['$\\sim$','$\\wedge$','n','s','Cfb','Dfb','Dfc'],dims=['tag'],
                        coords=[['shore','mountain','north','south','cfb','dfb','dfc']])
    print(f"########## Metrics Table {typ} #########")
    if not difdir:
        print(str("selection & N & "+
                  "$\\overline{\\mathrm{%s}}_\\mathrm{clr}$ & "%(typ)+
                  "$\\overline{\\mathrm{%s}}_\\mathrm{obs}$ & "%(typ)+
                  "$\\overline{\\mathrm{%s}}_\\mathrm{aer}$ & "%(typ)+
                  "RMSE & "+
                  "MBE & "+
                  "R &"+
                  "R$_\\mathrm{-clr}$ &"+
                  "$\\Delta\\!\\sigma_\\mathrm{all}$ "+
                  "\\\\"))
    if difdir:
        print(str("selection & N & "+
              "$\\overline{\\mathrm{%s}}_\\mathrm{clr}$ & "%(typ)+
              "$\\overline{\\mathrm{%s}}_\\mathrm{obs}$ & "%(typ)+
              "$\\overline{\\mathrm{%s}}_\\mathrm{aer}$ & "%(typ)+
              "$\\frac{\\overline{\\mathrm{DHI}}_\\mathrm{obs}}{\\overline{\\mathrm{DNI}\\mu_0}_\\mathrm{obs}}$ & "+
              "$\\frac{\\overline{\\mathrm{DHI}}_\\mathrm{aer}}{\\overline{\\mathrm{DNI}\\mu_0}_\\mathrm{aer}}$ & "+
              "RMSE & "+
              "MBE & "+
              "R &"+
              "R$_\\mathrm{-clr}$ &"+
              "$\\Delta\\!\\sigma_\\mathrm{all}$ "+
              "\\\\"))
        
    print('\\middlehline')
    for seas in ['MAM','JJA','SON','DJF']:        
        line = seas + ' &'
        line+=_add_data(dsy,typ=typ,st=None,seas=seas,tag=None,sza=sza,difdir=difdir)
        print(line)
    print('\\middlehline')
    for tag in ['shore','mountain','north','south','cfb','dfb',None]:
        if not tag ==None:
            if tag in ['north','south','shore','mountain']:
                line = tag + ' ({})'.format(str(tags.sel(tag=tag).values)) + ' &'
            else:
                line=str(tags.sel(tag=tag).values)+' &'
        else:
            print('\\middlehline')
            line='all &'
        line+=_add_data(dsy,typ=typ,st=None,seas=None,tag=tag,sza=sza,difdir=difdir)
        print(line)


    


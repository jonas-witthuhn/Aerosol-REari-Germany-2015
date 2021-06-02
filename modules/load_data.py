#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 13:10:02 2019

@author: walther
"""
import os
import json

import h5py
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt

import trosat.sunpos as sp
import modules.msevi_geolocation as msevi

def aeronet_stations():
    """ Store Metadata of all AERONET station in xarray Datset"""
    with open("datasets/aeronet_sites.json") as txt:
        sites=json.load(txt)
    return sites

def dwd_stations():
    """Store Metadata of all DWD stations in xarray Dataset"""
    stations={"AK":{'name':"Arkona",'coords':(54.68,13.44),'height': 42},
              "BG":{'name':"Braunschweig",'coords':(52.29,10.45),'height': 88},
              "BN":{'name':"Bremen (FWW)",'coords':(53.05,8.79),'height': 5},
              "CH":{'name':"Chemnitz",'coords':(50.81,12.89),'height': 357},
              "DN":{'name':"Dresden-Klotzsche",'coords':(51.13,13.77),'height': 222},
              "FB":{'name':"Fichtelberg",'coords':(50.43,12.96),'height': 1213},
              "FL":{'name':"Fürstenzell",'coords':(48.55,13.36),'height': 476},
              "GZ":{'name':"Görlitz",'coords':(51.16,14.96),'height': 238},
              "HF":{'name':"Hamburg-Fuhlsbüttel",'coords':(53.64,10.00),'height': 16},
              "HP":{'name':"Hohenpeißenberg",'coords':(47.80,11.02),'height': 977},
              "KS":{'name':"Konstanz",'coords':(47.68 ,9.19),'height': 443},
              "LG":{'name':"Lindenberg (RAO)",'coords':(52.21,14.12),'height': 98},
              "LZ":{'name':"Leipzig-Holzhausen",'coords':(51.32,12.41),'height': 148},
              "NB":{'name':"Nürnberg (FWW)",'coords':(49.50,11.08),'height': 312},
              "NY":{'name':"Norderney",'coords':(53.71 ,7.15),'height': 13},
              "PG":{'name':"St.Peter-Ording",'coords':(54.33 ,8.60),'height': 5},
              "PT":{'name':"Potsdam",'coords':(52.38,13.06 ),'height': 81},
              "RO":{'name':"Rostock-Warnemünde",'coords':(54.18,12.08 ),'height': 4},
              "SG":{'name':"Schleswig",'coords':(54.53,9.55 ),'height': 43},
              "SN":{'name':"Seehausen",'coords':(52.89 ,11.73),'height': 21},
              "SR":{'name':"Saarbrücken (FWW)",'coords':(49.22,7.11),'height': 320},
              "SY":{'name':"Stuttgart-Schnarrenberg",'coords':(48.83,9.20),'height': 311},
              "TR":{'name':"Trier",'coords':(49.75,6.66),'height': 265},
              "WN":{'name':"Weihenstephan",'coords':(48.40,11.73),'height': 467},
              "WZ":{'name':"Würzburg",'coords':(49.77,9.96),'height': 268},
              "ZG":{'name':"Zugspitze",'coords':(47.42,10.99),'height':2960 }}
    
    shoresite=['NY','PG','HF','SG','RO','AK']
    CFBsite = ['HF','BN','PG','SG','NY','RO','AK','SN','BG','TR','SR','SY']
    DFBsite = ['PT','LG','GZ','DN','CH','LZ','WZ','NB','FL','WN','HP','KS']
    DFCsite = ['ZG','FB']

    
    ## Mountain Site - Site is heigher than 400 m
    mountain_site=[stations[s]['height']>400. for s in stations.keys()]
    ## Shore Site - handselected
    shore_site=[True if s in shoresite else False for s in stations.keys()]
    #north_site=[True if s in northsite else False for s in stations.keys()]
    north_site=[stations[s]['coords'][0]>52 for s in stations.keys()]
    ## South site -  latitude < 50
    south_site=[stations[s]['coords'][0]<50 for s in stations.keys()]
    ## climate Cfb 
    cfb_site=[True if s in CFBsite else False for s in stations.keys()]
    ## climate Dfb 
    dfb_site=[True if s in DFBsite else False for s in stations.keys()]
    ## climate Dfc 
    dfc_site=[True if s in DFCsite else False for s in stations.keys()]
    
    
    tag=np.vstack((mountain_site,
		   shore_site,
		   north_site,
		   south_site,
		   cfb_site,
		   dfb_site,
		   dfc_site)).T
    
    
    ds=xr.Dataset({'name':('station',[stations[station]['name'] for station in stations.keys()]),
                   'latitude':('station',[stations[station]['coords'][0] for station in stations.keys()]),
                   'longitude':('station',[stations[station]['coords'][1] for station in stations.keys()]),
                   'height':('station',[stations[station]['height'] for station in stations.keys()]),
                   'site':(('station','tag'),tag)},
                  coords={'station':('station',list(stations.keys())),
                          'tag':('tag',['mountain','shore','north','south','cfb','dfb','dfc'])})
    return ds




def load_esd(dates):
    """ Store earth-sun-distance for input days in xarray Dataset """
    # unique days
    dates=np.unique(dates.astype('datetime64[D]'))
    # earth sun distance
    ESD=sp.earth_sun_distance(dates)
    ds=xr.Dataset({'esd':('time',ESD)},
                  coords={'time':('time',dates)})
    return ds


def load_tsi(pf="datasets/sorce_tsi_24hr_l3_full.csv"):
    """ Read total solar irradaince from SORCE and store in xarray Dataset. """
    df=pd.read_csv(pf,parse_dates=[0],date_parser=lambda x: x.astype("datetime64[D]"))
    ds=xr.Dataset.from_dataframe(df)
    for key in ds.keys():
        newkey=key.split()[0]
        ds=ds.rename_vars({key:newkey})
    ds=ds.swap_dims({'index':'time'})
    ds = ds.reindex({'time':pd.date_range("2003-01-01","2019-12-31")},method='nearest')
    return ds

def load_albedo(date,coords,pf):
    """
    Load albedo from lsasaf dataset    

    Parameters
    ----------
    date : np.datetime64
        single day to load the data from
    coords : np.ndarray(n,2), tuple (lat,lon), list [lat,lon]
        Coordinates to load the data.
        must be one of (float(lat),float(lon)),
                        (list(lat),list(lon)),
                        np.meshgrid(lat,lon),
                        np.array([[lat,lon],...,[lat,lon]])
    pf : string, optional
        path to lsasaf dataset of shape {pf}/{year}/{month}/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_{year}{month}{day}0000.
        The default is "/vols/satellite/home/jonas/workdata/lsasaf/".

    Returns
    -------
    albedoDH : np.ndarray of shape like lat and lon
        directional hemispheric albedo (black sky albedo)
    albedoBH : np.ndarray of shape like lat and lon
        bi-directional hemispheric albedo (white sky albedo)

    """
    
    ### make filepath string
    year,month,day=pd.to_datetime(date).strftime("%Y,%m,%d").split(',')
    name="{year}/{month}/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_{year}{month}{day}0000"
    fname=os.path.join(pf,name.format(year=year,month=month,day=day))
    if not os.path.exists(fname):
        return False,False
    ### translate coordinates
    if type(coords) in [list,tuple]:
        xx,yy=msevi.latlon2xy(np.array(coords[0]), np.array(coords[1]))
    elif type(coords) == np.ndarray:
        xx,yy=msevi.latlon2xy(coords[:,0],coords[:,1])
    else:
        raise TypeError("Incorrect type of 'coords' - must be one of (float(lat),float(lon)) or (list(lat),list(lon)) or np.meshgrid(lat,lon) or np.array([[lat,lon],...,[lat,lon]])")
    xx=np.round(xx).astype(int)
    yy=np.round(yy).astype(int)
    
    ### load data from daily file
    with h5py.File(fname,'r') as f:
        # white sky albedo - only diffuse light
        albedoBH=np.array(np.array(f['AL-BB-BH'])[yy,xx] / f['AL-BB-BH'].attrs['SCALING_FACTOR'])
        albedoBH[albedoBH==f['AL-BB-BH'].attrs['MISSING_VALUE']]=np.nan
        # black sky albedo - only direct light
        albedoDH=np.array(np.array(f['AL-BB-DH'])[yy,xx] / f['AL-BB-DH'].attrs['SCALING_FACTOR'])
        albedoDH[albedoDH==f['AL-BB-DH'].attrs['MISSING_VALUE']]=np.nan

    return albedoDH.T,albedoBH.T


def load_aeronet(site,pf):
    fname = "{site}.{sfx}"
    # load level 2.0 for aod and AE
    sfx_aod='{lvl}'.format(lvl='lev20')
    # load level 1.5 for inversion products
    # as otherwise dataset will be very sparse
    sfx_all='all.{lvl}'.format(lvl='lev15')

    # load aod
    try:
        f=os.path.join(pf,fname.format(site=site,sfx=sfx_aod))
        aer_aod = pd.read_csv(f,sep=',',header=6,parse_dates={'time':[0,1]},
                              date_parser=lambda x: pd.to_datetime(x,format='%d:%m:%Y %H:%M:%S'))
        aer_aod[aer_aod.values==-999]=np.nan
        aer_aod = aer_aod.dropna(how='all')
        aer_aod = aer_aod[aer_aod['time']>pd.to_datetime("2003-01-01")]
        if len(aer_aod)==0:
            aer_aod=False
        else:
            aer_aod = aer_aod.reset_index(drop=True)
    except:
        aer_aod=False
        
    # load almucantar product
    try:
        f=os.path.join(pf,fname.format(site=site,sfx=sfx_all))
        aer_all = pd.read_csv(f,sep=',',header=6,parse_dates={'time':[1,2]},
                              date_parser=lambda x: pd.to_datetime(x,format='%d:%m:%Y %H:%M:%S'))
        aer_all[aer_all.values==-999]=np.nan
        aer_all = aer_all.dropna(how='all')
        aer_all = aer_all[aer_all['time']>pd.to_datetime("2003-01-01")]
        if len(aer_all)==0:
            aer_all=False
        else:
            aer_all = aer_all.reset_index(drop=True)
    except:
        aer_all=False
    return aer_aod,aer_all


def closest_time(X,Y,dt=np.timedelta64(90,'m'),removedouble=True):
    """
    find best matching index of time within +- dt [s]
    Xdat - 1d array(N)
    Ydat - 1d array(M)
    dt - np.timedelta64
    
    output
    xind  - 1d array indexes (J)
    yind  - 1d array indexes (J)
    -> Xdat[xind] most equal to Ydat[yind]  (within dt limits)
    """
    switch=False
    if len(X)>len(Y):
        switch=True
    
    if switch:
        Xdat=Y.copy()
        Ydat=X.copy()
    else:
        Xdat=X.copy()
        Ydat=Y.copy()
        
    ids=np.searchsorted(Ydat,Xdat-dt) 
    xind=[]
    yind=[]
    
    
    if removedouble: # choose only uniqe pairs (closest match), eg. X[1] <->Y[10] and not X[1],X[2] <-> Y[10]
        for i in range(len(Xdat)):
            t1=Xdat[i]
            tempt=np.timedelta64(1,'Y').astype('timedelta64[ns]')
            new=True
            for j in range(len(Ydat[ids[i]:])):
                t2=Ydat[j+ids[i]]
                tmpdt=np.abs(t1-t2)
                if tmpdt<=dt and tmpdt< tempt :
                    if new:
                        xind.append(i)
                        yind.append(j+ids[i])
                        new=False
                    else:
                        xind[-1]=i
                        yind[-1]=j+ids[i]
                    tempt=tmpdt
                if tmpdt>tempt:
                    break
                if t2>t1+dt:
                    break
        xind,yind=np.array(xind),np.array(yind)
        i=0
        ilist=[]
        #remove doubles in Ydata
        while i < len(yind):
            tx,ty=Xdat[xind[yind==yind[i]]],Ydat[yind[i]]
            dtx=np.abs(np.array(tx,dtype=np.int)-int(ty))
            ilist.append(i+np.argmin(dtx))
            i+=len(dtx)
        xind,yind=xind[ilist],yind[ilist]
                
    else:
        for i in range(len(Xdat)):
            t1=Xdat[i]
            for j in range(len(Ydat[ids[i]:])):
                t2=Ydat[j+ids[i]]
                tmpdt=np.abs(t1-t2)
                if tmpdt<=dt:
                    #print i,t1,j,t2
                    xind.append(i)
                    yind.append(j+ids[i])
                if t2>t1+dt:
                    break  
    xind,yind=np.array(xind),np.array(yind)           
    if switch:
        return yind,xind
    else:
        return xind,yind
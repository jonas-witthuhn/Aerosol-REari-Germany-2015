#!/usr/bin/env python3
import os
import xarray as xr
import numpy as np
import configparser

import modules.load_data as ld


config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']
pf_lsasaf = os.path.join(pf,"LSASAF/")
dataset_albedo = os.path.join(pf,"ALBEDO.nc")

stations = ld.dwd_stations()
dates = np.arange("2015-01-01","2016-01-01",dtype='datetime64[D]')
lon=[]
lat=[]
for i,station in enumerate(stations.station.data):
    lat.append(stations.sel(station=station).latitude.data)
    lon.append(stations.sel(station=station).longitude.data)
lat=np.array(lat)
lon=np.array(lon)
for i,date in enumerate(dates):
    print(f"make albedo dataset for day: ",date)
    # load albedo
    ablack,awhite=ld.load_albedo(date, (lat,lon), pf_lsasaf)
    # skip day if no albedo data
    if type(ablack)==bool:
        ablack,awhite = np.zeros(len(lon))*np.nan,np.zeros(len(lon))*np.nan
    if i==0:
        ABLACK = np.array(ablack)
        AWHITE = np.array(awhite)
    else:
        ABLACK = np.vstack((ABLACK,np.array(ablack)))
        AWHITE = np.vstack((AWHITE,np.array(awhite)))    


ds = xr.Dataset({'ablack':(('time','station'),ABLACK),
                 'awhite':(('time','station'),AWHITE)},
                coords={'time': ('time',dates),
                        'station': ('station',stations.station.data)})

ds.to_netcdf(dataset_albedo,
             encoding=dict(ablack={'zlib':True},
                           awhite={'zlib':True}))
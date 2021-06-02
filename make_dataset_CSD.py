#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:27:33 2020

@author: walther
"""
import os
import xarray as xr
import numpy as np
import configparser

import modules.load_data as ld
import trosat.sunpos as sp
from clear_sky_models import models as csm
from clear_sky_detection.csd_bright_sun_2020 import bright_sun_csdc as bs_csdc
from clear_sky_detection.csd_bright_sun_2020 import bright_sun_csds as bs_csds


config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']

# path to DWD observations
pf_DWD = os.path.join(pf, "DWD/2015_{station}.nc")
# This dataset will be produced here, it stores the clear sky flags for 
# the DWD observations
dataset_CSD=os.path.join(pf,'CSD/BS_{station}.nc')

stations = ld.dwd_stations()

for station in stations.station.values:
    print(f"### {stations.sel(station=station).name.values}")
    data = xr.load_dataset(pf_DWD.format(station=station))
    zen,azi = sp.sun_angles(data.time.data,
                            data.latitude.data[0],
                            data.longitude.data[0])

    days = data.time.data.astype("datetime64[D]")
    
    ### calculate clear sky from model as a-priori input to clear sky detection 
    mdni,mdhi,mghi=[],[],[]
    for i,d in enumerate(np.unique(days)):
        ind=days==d
        mdni1,mdhi1,mghi1=csm.model_15_kasm(sza=data.solar_zenith[ind],
                                             date=d,
                                             pressure=data.surf_p_minute[ind],
                                             wv=data.tcwv_atmcm[ind])
        mdni.extend(mdni1)
        mdhi.extend(mdhi1)
        mghi.extend(mghi1)
        
    mdni = np.array(mdni)
    mdhi = np.array(mdhi)
    mghi = np.array(mghi)
    
    ### calculate clearsky detection
    ## csd  - cloud free
    csdc,_,_,bsdni,bsdhi = bs_csdc(data.time.data, data.solar_zenith.data,
                                   data['global'].data, mghi,
                                   data.diffuse_sc.data, mdhi,
                                   data.longitude.data[0])
    bsghi = bsdhi + bsdni*np.cos(np.deg2rad(data.solar_zenith.data))
    csdc[data.solar_zenith>90] = False
    csdc[csdc.mask] = False
    csdc = csdc.data
    ## csd - free-sun
    csds,_,_,_,_ = bs_csds(data.time.data, data.solar_zenith.data,
                           data['global'].data, mghi,
                           data.diffuse_sc.data, mdhi,
                           data.longitude.data[0])
    csds[data.solar_zenith>90] = False
    csds[csds.mask]=False
    csds = csds.data

    ## store results in xarray Dataset
    encoding = {'csdc': {'dtype':'bool','zlib':True},
                'csds': {'dtype':'bool','zlib':True}}


    dataset = xr.Dataset({'csdc' :("time",csdc),
                          'csds':("time",csds)},
                           coords={"time":('time',data.time.data)},
                           attrs={"clear_sky_model":"KASM",
                                  "clear_sky_detection":'Bright_Sun_2020'})
    dataset.to_netcdf(dataset_CSD.format(station=station),
                      encoding=encoding)


import os
import datetime as dt
import numpy as np
import xarray as xr
import pandas as pd
from scipy.optimize import minimize_scalar
import configparser

from trosat import sunpos as sp
import modules.load_data as ld
from clear_sky_models import models as csm
from modules import FourClassValidation 
FCV=FourClassValidation.FourClassValidation
normalize_metrics=FourClassValidation.normalize_metrics
os.nice(10)

# paths
config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']
pf_csd = os.path.join(pf,'CSD/BS_{station}.nc')
pf_csf = os.path.join(pf,'CSF/{model}_{station}.nc')
pf_dwd = os.path.join(pf,'DWD/2015_{station}.nc')
pf_albedo = os.path.join(pf,"ALBEDO.nc")
pf_tsi = os.path.join(pf,"sorce_tsi_24hr_l3_full.csv")


pf_save = os.path.join(pf,"test_data/CSFselection/{selection}_{station}_{model}_{dayofyear}.nc")
# iteration variables
days = pd.date_range('2015-01-01','2015-12-31')
stations= ld.dwd_stations()
models = ['MRM61',
          'MMAC',
          'Heliosat1I',
          'CEM',
          'ESRA',
          'METSTAT',
          'SOLISsimple']

def calc_albedo(R,awhite,ablack):
    return (1.-R)*ablack + R*awhite

def aggregate(csd):
    csdsum = csd.groupby('time.hour').sum()>1
    csdh=csdsum.astype(int)
    return csdh

def recursion_fkt(X,Xkey,GHI,CSM,kwargs):
    kwargs.update({Xkey:X})
    mdni,mdhi,mghi=CSM(**kwargs)
    vs=normalize_metrics(FCV(GHI,mghi,['A']))
    metricsum=0.
    for key in vs.dtype.names:
        if key in ['Om','Pm','N','TS','SBF','U95','R2']:
            continue
        metricsum+=np.nansum(vs[key])
    return metricsum

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
    # sum/mean
    SMO = np.sum(O)/np.mean(O)
    SMP = np.sum(P)/np.mean(O)
    # correlation
    R=np.corrcoef(O,P)[0,1]
    # MBE
    MBE = N1 * np.sum(DELTA)
    # RMSE
    RMSE = np.sqrt(N2*np.sum(DELTA**2))
    return N,SMO,SMP,R,MBE,RMSE


S0 = ld.load_tsi(pf_tsi)
### load albedo
ALBEDO = xr.load_dataset(pf_albedo)

encoding = dict(MBE={'zlib':True},
                RMSE={'zlib':True},
                mu0N={'zlib':True},
                mu0max={'zlib':True},
                mu0min={'zlib':True})

for station in stations.station.values:
    print('#############################################################')
    print('### ',station)
    CSD1 = xr.open_dataset(pf_csd.format(station=station))
    # choose longest clear sky day
#     dayofyear = np.argmax(CSD.csdc.groupby('time.dayofyear').sum().values)
    csdsum = CSD1.csdc.groupby('time.dayofyear').sum()
    dayofyears = csdsum.dayofyear.where(csdsum>300,drop=True).values
    for dayofyear in dayofyears:
        print('#############################################################')
        print('### ',dayofyear)
        day = dt.datetime(2015,1,1)+dt.timedelta(int(dayofyear))
        DWD = xr.open_dataset(pf_dwd.format(station=station))
        CSD = xr.open_dataset(pf_csd.format(station=station))
        DWD = DWD.sel(time=day.strftime("%Y-%m-%d"))
        CSD = CSD.sel(time=day.strftime("%Y-%m-%d"))

        # calculate albedo 
        tmpALBEDO = ALBEDO.interp(time=DWD.time.values,method='linear')
        ablack = tmpALBEDO.sel(station = station).ablack
        awhite = tmpALBEDO.sel(station = station).awhite
        #diffuse to direct ratio
        Rdwd=DWD.diffuse_sc/(DWD['global']-DWD.diffuse_sc) 
        #albedo
        albedo = calc_albedo(Rdwd,awhite,ablack)


        DWD = DWD.where((0<=Rdwd)*(Rdwd<=1),drop=True)
        CSD = CSD.where((0<=Rdwd)*(Rdwd<=1),drop=True)
        albedo = albedo.where((0<=Rdwd)*(Rdwd<=1),drop=True)

        # get hour where sky is clear more than 2min
        csdhours = CSD.csdc.groupby('time.hour').sum()>2
        HoursOfSun = csdhours.hour.where(csdhours,drop=True).values

    #     for model in models:
        for model in ['ESRA']:
            print('#################################################')
            print('### ',model)

            MBE1h=[]
            RMSE1h=[]
            mu0N1h=[]
            mu0max1h=[]
            mu0min1h=[]
            print("---- single hours -------")
            for hour in HoursOfSun:
                idx = (CSD.csdc.where(CSD.time.dt.hour == hour,other=False)).astype(bool)
                dwd = DWD.where(idx,drop=True)
                mu0range = (np.max(dwd.solar_zenith.values),np.min(dwd.solar_zenith.values))
                mu0N = len(np.unique(np.round(dwd.solar_zenith.values,0)))
                # daily mean albedo
                A = np.mean(albedo.where(idx,drop=True).values)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)

                GHI=dwd['global'].data
                args=('aod',
                      GHI,
                      csm.clear_sky_model_with_aod,
                      csm_kwargs.copy())
                res=minimize_scalar(recursion_fkt,bounds=(0.,1.5),args=args)

                ## calculate optimized clearsky for whole day
                dwd = DWD.where(CSD.csdc,drop=True)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)
                time=dwd.time.data
                aod0=np.ones(len(time))*res.x
                # calculate with optimal aod
                mdni,mdhi,mghi=csm.clear_sky_model_with_aod(aod = res.x, **csm_kwargs)
                N,SMO,SMP,R,MBE,RMSE = get_metrics(dwd['global'].values,mghi)
                MBE1h.append(MBE)
                RMSE1h.append(RMSE)
                mu0N1h.append(mu0N)
                mu0max1h.append(mu0range[1])
                mu0min1h.append(mu0range[0])
                print(f"Hour ({hour:.0f}) - szen range: {mu0N} in ({mu0range[0]:.0f},{mu0range[1]:.0f}) - {R:.2f}, {MBE: 7.2f}, {RMSE:.2f}")
            # save     
            ds = xr.Dataset(dict(MBE = np.array(MBE1h),
                                   RMSE = np.array(RMSE1h),
                                   mu0N = np.array(mu0N1h),
                                   mu0max = np.array(mu0max1h),
                                   mu0min = np.array(mu0min1h)))
            ds.to_netcdf(pf_save.format(selection='1h',station=station,model=model,dayofyear=dayofyear),
                         encoding=encoding)

            print("---- two hours -------")
            # no repeated combination of two of clear sky hours
            newhours = np.array(np.meshgrid(HoursOfSun,HoursOfSun)).T.reshape(-1, 2)
            HoursOfSun2 = np.array([list(a[:])  for a in newhours if (a[0]<a[1])])
            MBE2h=[]
            RMSE2h=[]
            mu0N2h=[]
            mu0max2h=[]
            mu0min2h=[]
            for hours in HoursOfSun2:
                selection=(CSD.time.dt.hour==hours[0])+(CSD.time.dt.hour == hours[1])
                idx = (CSD.csdc.where(selection,other=False)).astype(bool)
                dwd = DWD.where(idx,drop=True)
                mu0range = (np.max(dwd.solar_zenith.values),np.min(dwd.solar_zenith.values))
                mu0N = len(np.unique(np.round(dwd.solar_zenith.values,0)))
                # daily mean albedo
                A = np.mean(albedo.where(idx,drop=True).values)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)

                GHI=dwd['global'].data
                args=('aod',
                      GHI,
                      csm.clear_sky_model_with_aod,
                      csm_kwargs.copy())
                res=minimize_scalar(recursion_fkt,bounds=(0.,1.5),args=args)

                ## calculate optimized clearsky for whole day
                dwd = DWD.where(CSD.csdc,drop=True)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)
                time=dwd.time.data
                aod0=np.ones(len(time))*res.x
                # calculate with optimal aod
                mdni,mdhi,mghi=csm.clear_sky_model_with_aod(aod = res.x, **csm_kwargs)
                N,SMO,SMP,R,MBE,RMSE = get_metrics(dwd['global'].values,mghi)
                MBE2h.append(MBE)
                RMSE2h.append(RMSE)
                mu0N2h.append(mu0N)
                mu0max2h.append(mu0range[1])
                mu0min2h.append(mu0range[0])
                print(f"Hours ({hours[0]:.0f},{hours[1]:.0f}) - szen range: {mu0N} in ({mu0range[0]:.0f},{mu0range[1]:.0f}) - {R:.2f}, {MBE: 7.2f}, {RMSE:.2f}")

            # save     
            ds = xr.Dataset(dict(MBE = np.array(MBE2h),
                                   RMSE = np.array(RMSE2h),
                                   mu0N = np.array(mu0N2h),
                                   mu0max = np.array(mu0max2h),
                                   mu0min = np.array(mu0min2h)))
            ds.to_netcdf(pf_save.format(selection='2h',station=station,model=model,dayofyear=dayofyear),
                         encoding=encoding)

            print("---- three hours -------")
            # no repeated combination of two of clear sky hours
            newhours = np.array(np.meshgrid(HoursOfSun,HoursOfSun,HoursOfSun)).T.reshape(-1, 3)
            HoursOfSun3 = np.array([list(a[:])  for a in newhours if (a[0]<a[1]<a[2])])
            MBE3h=[]
            RMSE3h=[]
            mu0N3h=[]
            mu0max3h=[]
            mu0min3h=[]
            for hours in HoursOfSun3:
                selection=(CSD.time.dt.hour==hours[0])
                selection+=(CSD.time.dt.hour == hours[1])
                selection+=(CSD.time.dt.hour == hours[2])
                idx = (CSD.csdc.where(selection,other=False)).astype(bool)
                dwd = DWD.where(idx,drop=True)
                mu0range = (np.max(dwd.solar_zenith.values),np.min(dwd.solar_zenith.values))
                mu0N = len(np.unique(np.round(dwd.solar_zenith.values,0)))
                # daily mean albedo
                A = np.mean(albedo.where(idx,drop=True).values)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)

                GHI=dwd['global'].data
                args=('aod',
                      GHI,
                      csm.clear_sky_model_with_aod,
                      csm_kwargs.copy())
                res=minimize_scalar(recursion_fkt,bounds=(0.,1.5),args=args)

                ## calculate optimized clearsky for whole day
                dwd = DWD.where(CSD.csdc,drop=True)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)
                time=dwd.time.data
                aod0=np.ones(len(time))*res.x
                # calculate with optimal aod
                mdni,mdhi,mghi=csm.clear_sky_model_with_aod(aod = res.x, **csm_kwargs)
                N,SMO,SMP,R,MBE,RMSE = get_metrics(dwd['global'].values,mghi)
                MBE3h.append(MBE)
                RMSE3h.append(RMSE)
                mu0N3h.append(mu0N)
                mu0max3h.append(mu0range[1])
                mu0min3h.append(mu0range[0])
                print(f"Hours ({hours[0]:.0f},{hours[1]:.0f},{hours[2]:.0f}) - szen range: {mu0N} in ({mu0range[0]:.0f},{mu0range[1]:.0f}) - {R:.2f}, {MBE: 7.2f}, {RMSE:.2f}")

            # save     
            ds = xr.Dataset(dict(MBE = np.array(MBE3h),
                                   RMSE = np.array(RMSE3h),
                                   mu0N = np.array(mu0N3h),
                                   mu0max = np.array(mu0max3h),
                                   mu0min = np.array(mu0min3h)))
            ds.to_netcdf(pf_save.format(selection='3h',station=station,model=model,dayofyear=dayofyear),
                         encoding=encoding)            

            print("---- four hours -------")
            # no repeated combination of two of clear sky hours
            newhours = np.array(np.meshgrid(HoursOfSun,HoursOfSun,HoursOfSun,HoursOfSun)).T.reshape(-1, 4)
            HoursOfSun4 = np.array([list(a[:])  for a in newhours if (a[0]<a[1]<a[2]<a[3])])
            MBE4h=[]
            RMSE4h=[]
            mu0N4h=[]
            mu0max4h=[]
            mu0min4h=[]
            for hours in HoursOfSun4:
                selection=(CSD.time.dt.hour==hours[0])
                selection+=(CSD.time.dt.hour == hours[1])
                selection+=(CSD.time.dt.hour == hours[2])
                selection+=(CSD.time.dt.hour == hours[3])
                idx = (CSD.csdc.where(selection,other=False)).astype(bool)
                dwd = DWD.where(idx,drop=True)
                mu0range = (np.max(dwd.solar_zenith.values),np.min(dwd.solar_zenith.values))
                mu0N = len(np.unique(np.round(dwd.solar_zenith.values,0)))
                # daily mean albedo
                A = np.mean(albedo.where(idx,drop=True).values)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)

                GHI=dwd['global'].data
                args=('aod',
                      GHI,
                      csm.clear_sky_model_with_aod,
                      csm_kwargs.copy())
                res=minimize_scalar(recursion_fkt,bounds=(0.,1.5),args=args)

                ## calculate optimized clearsky for whole day
                dwd = DWD.where(CSD.csdc,drop=True)
                csm_kwargs=dict(name=model,
                                sza=dwd.solar_zenith.values,
                                date=np.datetime64(day).astype('datetime64[D]'),
                                pressure=dwd.surf_p_minute.values,
                                altitude = stations.sel(station=station).height.values,
                                wv=dwd.tcwv_atmcm.values,
                                ozone=dwd.tco3_atmcm.values,
                                albedo=A,
                                Esc=S0.sel(time=day.strftime('%Y%m%d')).tsi_1au.values)
                time=dwd.time.data
                aod0=np.ones(len(time))*res.x
                # calculate with optimal aod
                mdni,mdhi,mghi=csm.clear_sky_model_with_aod(aod = res.x, **csm_kwargs)
                N,SMO,SMP,R,MBE,RMSE = get_metrics(dwd['global'].values,mghi)
                MBE4h.append(MBE)
                RMSE4h.append(RMSE)
                mu0N4h.append(mu0N)
                mu0max4h.append(mu0range[1])
                mu0min4h.append(mu0range[0])
                print(f"Hours ({hours[0]:.0f},{hours[1]:.0f},{hours[2]:.0f},{hours[3]:.0f}) - szen range: {mu0N} in ({mu0range[0]:.0f},{mu0range[1]:.0f}) - {R:.2f}, {MBE: 7.2f}, {RMSE:.2f}")

            # save     
            ds = xr.Dataset(dict(MBE = np.array(MBE4h),
                                   RMSE = np.array(RMSE4h),
                                   mu0N = np.array(mu0N4h),
                                   mu0max = np.array(mu0max4h),
                                   mu0min = np.array(mu0min4h)))
            ds.to_netcdf(pf_save.format(selection='4h',station=station,model=model,dayofyear=dayofyear),
                         encoding=encoding)


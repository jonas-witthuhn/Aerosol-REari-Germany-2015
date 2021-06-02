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

S0 = ld.load_tsi(pf_tsi)

# lookup timeindex
ds = xr.load_dataset(pf_dwd.format(station='LG'))
TIME = ds.time.data
del ds

### load albedo
ALBEDO = xr.load_dataset(pf_albedo)
ALBEDO = ALBEDO.interp(time=TIME,method='linear')

AREdata_dwd = np.zeros((len(days),len(stations.station),len(models)))
AREdata_csf = np.zeros((len(days),len(stations.station),len(models)))

for m,model in enumerate(models):
    for s,station in enumerate(stations.station.data):
        ds_flag=False
        print(f'CSF with {model} ({m+1}/{len(models)}) at {station} ({s+1}/{len(stations.station.data)})')
        DWD = xr.load_dataset(pf_dwd.format(station=station))
        CSD = xr.load_dataset(pf_csd.format(station=station))

        # calculate albedo 
        ablack = ALBEDO.sel(station = station).ablack
        awhite = ALBEDO.sel(station = station).awhite
        #diffuse to direct ratio
        Rdwd=DWD.diffuse_sc/(DWD['global']-DWD.diffuse_sc) 
        #albedo
        albedo = calc_albedo(Rdwd,awhite,ablack)

        # select only if 0<R<1 and cloudfree for observation
        idx = CSD.csds.values*(0<=Rdwd)*(Rdwd<=1)

        # consecutive clearsky hours of all days
        csd = CSD.csdc
        csdh = csd.groupby('time.dayofyear').map(aggregate)
        csdh = np.sum(csdh,axis=1)

        # return days for potential clear sky fit
        # and counts of clear sky minutes identified during 
        # this day
        days,counts = np.unique(DWD.where(idx,drop=True).time.astype('datetime64[D]'),
                                return_counts=True)
        
        # CSF requirement:
        #  - at least 10 min of clear sky at the day
        #  - clear sky in two different hours of the day
        days = days[counts>10]
        days = days[csdh[pd.to_datetime(days).dayofyear-1] > 2]
        days = pd.to_datetime(days)
        for d,day in enumerate(days):
            print(f"  >> CSF at {day:%d.%m.%Y} ({d+1}/{len(days)})")
            # apply selection
            dwd = DWD.where(idx,drop=True).sel(time=day.strftime('%Y%m%d'))
            # daily mean albedo
            A = np.mean(albedo.where(idx,drop=True).sel(time=day.strftime('%Y%m%d')).values)
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
            dwd = DWD.sel(time=day.strftime('%Y%m%d'))
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
            # calculate with "clean air"
            mdni0,mdhi0,mghi0=csm.clear_sky_model_with_aod(aod = 0., **csm_kwargs) 


            ## prepare xarray dataset
            ds=xr.Dataset({'aod0': ('time',aod0),
                             'mdni0':('time',mdni0),
                             'mdhi0':('time',mdhi0),
                             'mghi0':('time',mghi0),
                             'mdni':('time', mdni),
                             'mdhi':('time',mdhi),
                             'mghi':('time',mghi)},
                            coords={'time':('time',time)},
                            attrs={'clear_sky_model':model})
            
            if not ds_flag:
                ## initialize and save first time
                encoding={'mdni':{'zlib':True},
                          'mghi':{'zlib':True},
                          'mdhi':{'zlib':True},
                          'mdni0':{'zlib':True},
                          'mghi0':{'zlib':True},
                          'mdhi0':{'zlib':True},
                          'aod0': {'zlib':True}}
                mfdata=ds.copy()
                ds_flag=True
            else:
                ## merge
                mfdata= xr.combine_by_coords((mfdata,ds))
                
        ## save 
        mfdata.to_netcdf(path=pf_csf.format(model=model,station=station),
                         encoding=encoding) 

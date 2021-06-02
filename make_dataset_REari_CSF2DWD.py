import os
import numpy as np
import xarray as xr
import pandas as pd
from scipy.integrate import quad,trapz
import configparser

import modules.load_data as ld
os.nice(10)



# paths
config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']



# define paths to
# clear sky fit datasets
pf_csf = os.path.join(pf,'CSF/{model}_{station}.nc')
# clear sky detection datasets
pf_csd = os.path.join(pf,'CSD/BS_{station}.nc')
# DWD irradiance observations
pf_dwd = os.path.join(pf,'DWD/2015_{station}.nc')
# LSA-SAF surface albedo 
pf_albedo = os.path.join(pf,"ALBEDO.nc")
# this dataset is produced here
dataset_REari_DWD = os.path.join(pf,"REari_CSF2DWD.nc")

# iteration variables
days = pd.date_range('2015-01-01','2015-12-31')
models = ['MRM61',
          'MMAC',
          'Heliosat1I',
          'CEM',
          'ESRA',
          'METSTAT',
          'SOLISsimple',
          ]
stations= ld.dwd_stations()


################################################################################################
### DWD + CSM
################################################################################################

## functions of calculating surface albeo and REari (ARE)
def calc_albedo(R,awhite,ablack):
    return (1.-R)*ablack + R*awhite
def calc_ARE(A,Faer,Fclr):
    return (1.-A)*(Faer-Fclr)


# lookup timeindex
ds = xr.load_dataset(pf_dwd.format(station='LG'))
TIME = ds.time.data
del ds

### load albedo and interpolate to DWD observations
ALBEDO = xr.load_dataset(pf_albedo)
ALBEDO = ALBEDO.interp(time=TIME,method='linear')


# initialize result arrays
AREdata_csf = np.zeros((len(days),len(stations.station),len(models)))
Adata_csf = np.zeros((len(days),len(stations.station),len(models)))
AREdata_csf_noal = np.zeros((len(days),len(stations.station),len(models)))


## calculate REari for each model+ DWD station combination
for m,model in enumerate(models):
    for s,station in enumerate(stations.station.data):
        if station == 'ZG':
            continue
        print(f'calculate REari for {model} ({m+1}/{len(models)}) - {station} ({s+1}/{len(stations.station.data)})')
        DWD = xr.load_dataset(pf_dwd.format(station=station))
        CSD = xr.load_dataset(pf_csd.format(station=station))
        CSF = xr.load_dataset(pf_csf.format(station=station,model=model))
        
        ## filter fit artifacts:
        ## Force 0 < AOD < 0.7, as due to CSF AOD sometimes jumps to negative or very large values
        CSF = CSF.where(CSF.aod0>0,drop=True)
        CSF = CSF.where(CSF.aod0<0.7,drop=True)
        ## clear sky fitting the MRM61 model leads to artifacts with values of exact 0.3, due to the
        ## selection utilizing aod in the model itself, therefore values with values of exact 0.3 are
        ## filtered out
        if model == 'MRM61':
            CSF = CSF.where(CSF.aod0.round(3)!=0.3,drop=True)
        
        # unify time index
        CSF = CSF.reindex_like(DWD)

        ### force fitted data to 0 if sun is not visible
        CSF = CSF.where(DWD.solar_zenith<=90,other=0) # set irradiance to 0
        CSF['aod0'] = CSF.aod0.where(DWD.solar_zenith<=90) # set aod to nan instead
        
        # calculate diffuse to direkt ratio for albedo calculation
        mu0 = np.cos(np.deg2rad(DWD.solar_zenith))
        Rcsf=CSF.mdhi/(CSF.mdni*mu0)
        
        # calculate albedo and REari
        ablack = ALBEDO.sel(station = station).ablack
        awhite = ALBEDO.sel(station = station).awhite
        Acsf = calc_albedo(Rcsf,awhite,ablack)
        AREcsf = calc_ARE(Acsf,CSF.mghi,CSF.mghi0)
        AREcsf_noal = calc_ARE(0.,CSF.mghi,CSF.mghi0)

        # index selector: select only if 0<R<1 and cloudfree for observation
        idxcsf = (0<=Rcsf)*(Rcsf<=1)

        for d,day in enumerate(days):
            time = day.strftime("%Y-%m-%d")
            # CSF daily averages
            X = CSF.sel(time=time).time[idxcsf.sel(time=time)]
            X = (X.dt.hour + X.dt.minute/60.)/24.
            try:
                VAL = AREcsf.where(idxcsf,drop=True).sel(time=time)
                AREdata_csf[d,s,m]=trapz(VAL,X)
            except:
                AREdata_csf[d,s,m]=np.nan
            
            try:
                VAL = AREcsf_noal.where(idxcsf,drop=True).sel(time=time)
                AREdata_csf_noal[d,s,m]=trapz(VAL,X)
            except:
                AREdata_csf_noal[d,s,m]=np.nan
            
            try:
                VAL = Acsf.where(idxcsf,drop=True).sel(time=time)
                Adata_csf[d,s,m]=trapz(VAL,X)
            except:
                Adata_csf[d,s,m]=np.nan
            ## optionally print daily REari 
#             if ~np.isnan(AREdata_csf[d,s,m]) and AREdata_csf[d,s,m]<0:
#                 print(time,np.round(AREdata_csf[d,s,m],2))
            
        del DWD,CSF,CSD
        del AREcsf,Acsf,idxcsf,Rcsf
        del mu0,ablack,awhite

Adata_csf[Adata_csf>=0]=np.nan
AREdata_csf[AREdata_csf>=0]=np.nan
AREdata_csf_noal[AREdata_csf_noal>=0]=np.nan

    

ds = xr.Dataset({'ARE_csf':(('day','station','model'),AREdata_csf),
                 'ARE_csf_noal':(('day','station','model'),AREdata_csf_noal),
                 'A_csf':(('day','station','model'),Adata_csf)},
                coords = {'day':('day',days),
                          'station':('station',stations.station.data),
                          'model':('model',models)})
encoding ={}
for key in ds.keys():
    encoding.update({key:dict(zlib=True)})
ds.to_netcdf(dataset_REari_DWD,
             encoding=encoding)
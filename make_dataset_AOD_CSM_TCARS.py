import os
os.nice(10)
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import configparser


import modules.load_data as ld
from CAMS_aerosol_props.load_cams_data import CAMS

#######################################################################
# paths
config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']

# The camsra aerosol mass mixing ratios will be used
# to calculate the spectral AOD at 550 and 700nm
# to compare to CSM
fname_camsra =  os.path.join(pf,"CAMSRA/cams-ra_{date:%Y-%m-%d}_{levtype}.nc")

# This dataset is produced by "make_dataset_ecrad_skill.py"
# In this dataset the simulated irradiance with and without aerosols
# by T-CARS is stored, collocated to the DWD stations.
# From this, the T-CARS broadband AOD will be calculated
fname_TCARS_irradiance = os.path.join(pf,"ecrad_dwd_skill.nc")

# This dataset is produced by "make_dataset_CSF.py"
# The dataset stores the clear sky fitted irradiance with one CSM
# and the AOD inverted from the CSM
fname_CSF = os.path.join(pf,'CSF/{model}_{station}.nc')

# This dataset is produced here,
# it stores AOD at 550 and 700nm collocated to DWD stations
# calculated from CAMS RA mass mixing ratios
fname_CAMS_AOD = os.path.join(pf,"TCARS_CSF_AOD_2015.nc")


# The DWD irradiance observations
fname_DWD_irradiance = os.path.join(pf,'DWD/2015_{station}.nc')

##########################################################################
# settings
dates = pd.date_range("2015-01-01","2015-12-31")
stations = ld.dwd_stations()
# we drop Zugspitze, since data is not available complete 2015
stations = stations.drop_sel(station='ZG')
# CSM names
models = np.array(['MRM61',
                  'MMAC',
                  'Heliosat1I',
                  'CEM',
                  'ESRA',
                  'METSTAT',
                  'SOLISsimple'])


# initialyze result dataset
shape1 = (len(dates),len(stations.station.values))
shape2 = (len(dates),len(stations.station.values),len(models))
dummy1 = np.zeros(shape1)*np.nan
dummy2 = np.zeros(shape2)*np.nan
dsAOD = xr.Dataset({'AOD0':(('day','station'),dummy1.copy()),
                    'AOD550':(('day','station'),dummy1.copy()),
                    'AOD700':(('day','station'),dummy1.copy()),
                    'CAOD':(('day','station','model'),dummy2.copy())},
                   coords={'day':('day',dates),
                           'station':('station',stations.station.values),
                           'model':('model',models)})


##########################################################################
### add broadband AOD
for s,station in enumerate(stations.station.values):
    BB = xr.open_dataset(fname_TCARS_irradiance)
    BBsel = BB.where(BB.mu0>0,drop=True)
    BBsel = BBsel.where(BBsel.dni_clr>0,drop=True)
    BBsel = BBsel.where(BBsel.dni_aer>0,drop=True)
    BBsel = BBsel.where(BBsel.name==station,drop=True)
    BBAOD = -BBsel.mu0*np.log(BBsel.dni_aer/BBsel.dni_clr)
    BBAOD = BBAOD.groupby('time.dayofyear').mean(skipna=True)
    BBAOD = BBAOD.swap_dims({'dayofyear':'day'}).assign_coords({'day':dates[BBAOD.dayofyear-1]})
    BBAOD = BBAOD.reindex_like(dsAOD.sel(station=station).AOD550)
    dsAOD.AOD0.values[:,s] = BBAOD.values

    #################################################################
    ### add AOD from CSF
    for m,model in enumerate(models):
        CSF = xr.open_dataset(fname_CSF.format(model=model,station=station))
        CAOD = CSF.aod0.groupby('time.dayofyear').mean(skipna=True)
        CAOD = CAOD.swap_dims({'dayofyear':'day'}).assign_coords({'day':dates[CAOD.dayofyear-1]})
        CAOD = CAOD.reindex_like(dsAOD.sel(station=station).AOD550)
        # Filter artifacts from the clear sky fit
        CAOD = CAOD.where(CAOD>0)
        CAOD = CAOD.where(CAOD<0.7)
        CAOD = CAOD.where(CAOD.round(3)!=0.3)
        # store to array
        dsAOD.CAOD.values[:,s,m] = CAOD.values


####################################################################################
### Calculate spectral AOD at 550 and 700nm from CAMS data
for d,day in enumerate(dates):
    # skip if clear sky fittet aod is not available
    if np.all(np.isnan(dsAOD.CAOD.sel(day=day).values)):
        continue
    print(day)
    ## To scale CAMS gridded data to stations, the station altitude and measured pressure
    ## is taken from the DWD dataset
    ## The arrays have to be the shape of the CAMS data after lat/lon interpolation
    with xr.open_dataset(fname_camsra.format(date=day,levtype='sfc')) as cams_sfc:
        # get altitude for each station, repeated for each timestep
        dwdalt = np.repeat(stations.height.values[np.newaxis,:],
                           cams_sfc.time.size,axis=0)
        # get surface pressure at each station
        for i,station in enumerate(stations.station.values):
            # load DWD station data and select CAMS time
            with xr.open_dataset(fname_DWD_irradiance.format(station=station)) as dwd:
                dwdp0_tmp = dwd.sel(time = cams_sfc.time).surf_p_minute.values
                dwdp0_tmp*=1e2 #convert from [hPa] to [Pa]
            if i==0:
                dwdp0 = dwdp0_tmp.copy()
            else:
                dwdp0 = np.vstack((dwdp0,dwdp0_tmp))
        dwdp0 = dwdp0.T

    dwdlat = stations.latitude.values
    dwdlon = stations.longitude.values
    ## load CAMS data and interpolate/scale to DWD stations
    CamsData = CAMS(fname_camsra.format(date=day,levtype='sfc'),
                    fname_camsra.format(date=day,levtype='ml'),
                    interp_lat=dwdlat,
                    interp_lon=dwdlon,
                    scale_altitude=dwdalt,
                    scale_sfcpressure=dwdp0,
                    combine_latlon=True)

    ## calculate aerosol optical properties
    AerosolProps,_ = CamsData.aerosol_optprop([550,700],mono=False)

    ## store daily mean AOD to dataset
    newshape = (CamsData.cams_sfc.time.size,stations.station.size)
    AOD550 = AerosolProps.aod.sel(wvl=550).values.reshape(newshape)
    AOD550 = np.mean(AOD550,axis=0)
    AOD700 = AerosolProps.aod.sel(wvl=700).values.reshape(newshape)
    AOD700 = np.mean(AOD700,axis=0)
    
    dsAOD.AOD550.values[d,:] = AOD550
    dsAOD.AOD700.values[d,:] = AOD700

################################################################################
# save dataset
dsAOD.to_netcdf(fname_CAMS_AOD,
                encoding=dict(AOD550={'zlib':True},
                              AOD700={'zlib':True},
                              AOD0={'zlib':True},
                              CAOD={'zlib':True}))
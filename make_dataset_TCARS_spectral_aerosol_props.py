import os
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import configparser
from scipy.special import roots_legendre

# custom modules
import modules.load_data as ld
import ecrad_python.camsra2ecrad as c2e
import ecrad_python.run_ecrad as ecrad
import ecrad_python.AE_perturbation as perturbeAE
from CAMS_aerosol_props.load_cams_data import CAMS
from trosat import sunpos as sp

# set nice value, since calculation might take a while, and uses cpu alot
os.nice(10)

# paths
config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']
ncfile_path = os.path.join(pf,"CAMSRA/")

# aerosol optical properties of cams aerosol
# https://doi.org/10.24380/jgs8-sc58
# have to be stored in <ecrad_path>/data/
aerosol_props = "aerosol_cams_ifs_optics.nc"

# filename templates
camsra_fname = "cams-ra_{date:%Y-%m-%d}_{levtype}.nc"

# dataset of the total solar irradiance from SORCE
tsi_path = os.path.join(pf,"sorce_tsi_24hr_l3_full.csv")

# REari kernel, REari, surface irradiance and aerosol optical props
# will be stored in this file
fname_dsTCARS = os.path.join(pf,'test_TCARS_spectral_aerosol_props_{year}.nc')


years = np.arange(2003,2020)
for year in years:
#     if year == 2015:
#         continue
    # define timeframe
    days = pd.date_range(f"{year}-01-01",f"{year}-12-31")
    for d,day in enumerate(days):
        print(f"  >> {day:%d.%m.%Y}")
        ##########################################################################
        # make ecRad input file

        fn_cams_sfc = os.path.join(ncfile_path,
                                   camsra_fname.format(date=day,levtype='sfc'))
        fn_cams_ml = os.path.join(ncfile_path,
                                  camsra_fname.format(date=day,levtype='ml'))

        CAMS_sfc = xr.open_dataset(fn_cams_sfc)

        # reduce grid
        lats=CAMS_sfc.lat.values
        lons=CAMS_sfc.lon.values
        sel_lat = np.arange(np.min(lats),np.max(lats)+0.75,0.75)
        sel_lon = np.arange(np.min(lons),np.max(lons)+0.75,0.75)


        #######################################################################
        # load CAMS RA data and make ecRad input
        # interpolate to reduced grid and time selection
        CamsData = CAMS(fn_cams_sfc,
                        fn_cams_ml,
                        interp_lat=sel_lat,
                        interp_lon=sel_lon)

        cams_sfc = CamsData.cams_sfc
        cams_ml = CamsData.cams_ml

        ####################################################################
        # calculate bulk aerosol optical properties
        # at ecmwf mean bands
        AERCFG = xr.open_dataset("CAMS_aerosol_props/aerosol_cams_ifs_optics.nc")
        channels1 = 1./AERCFG.wavenumber1_sw[:-1]
        channels2 = 1./AERCFG.wavenumber2_sw[:-1]
        wvls_bands=np.mean(np.vstack((channels1,channels2)),axis=0)
        wvls_bands*= 1e7 #[nm]
        AerosolProps,_ = CamsData.aerosol_optprop(wvls_bands,mono=False,delta_eddington=True)

        # at mono wavelengths
        wvls_mono = [440,469,550,670,865,1240]
        AerosolPropsMono,_ = CamsData.aerosol_optprop(wvls_mono,mono=True,delta_eddington=True)

        APtmp = xr.merge((AerosolProps,AerosolPropsMono))

        ## calculate aerosol  column mass
        mass = np.zeros(cams_sfc.psfc.shape)
        for key in cams_sfc.keys():
            if key[:6] == 'aermss':
                mass+= cams_sfc[key].values * 1000.  # [g m-2]


        ####################################################################
        ## reshape to dataset
        times,ind = np.unique(APtmp.time,return_index=True)
        times = times[np.argsort(ind)]
        lats,ind = np.unique(APtmp.lat,return_index=True)
        lats = lats[np.argsort(ind)]
        lons,ind = np.unique(APtmp.lon,return_index=True)
        lons = lons[np.argsort(ind)]

        newshape=(len(times),
                  len(lats),
                  len(lons),
                  APtmp.wvl.size)
        coords = dict(time=('time',times),
                      lat=('lat',lats),
                      lon=('lon',lons),
                      wvl=('wvl',APtmp.wvl))

        ap = xr.Dataset(dict(aod = (('time','lat','lon','wvl'),
                                APtmp.aod.values.reshape(newshape)),
                             ssa = (('time','lat','lon','wvl'),
                                APtmp.ssa.values.reshape(newshape)),
                             g = (('time','lat','lon','wvl'),
                                APtmp.g.values.reshape(newshape)),
                             mass = (('time','lat','lon'),mass)
                            ),
                         coords = coords)

        if d==0:
            AP = ap.copy()
        else:
            AP = AP.merge(ap,compat='no_conflicts')


    AP.to_netcdf(fname_dsTCARS.format(year=year),
                 encoding=dict(time = {'zlib':True},
                                 lat = {'zlib':True},
                                 lon = {'zlib':True},
                                 g = {'zlib':True},
                                 ssa = {'zlib':True},
                                 aod = {'zlib':True},
                                 mass = {'zlib':True}))    

import os
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import configparser
import json
from scipy.integrate import trapz

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
lsasaf_albedo_path = os.path.join(pf,"LSASAF/")
# path to ecrad installation
# https://github.com/ecmwf/ecrad
ecrad_path = config['PATHS']['ecrad']

# aerosol optical properties of cams aerosol
# https://doi.org/10.24380/jgs8-sc58
# have to be stored in <ecrad_path>/data/
aerosol_props = "aerosol_cams_ifs_optics.nc"

# filename templates
camsra_fname = "cams-ra_{date:%Y-%m-%d}_{levtype}.nc"

# dataset of the total solar irradiance from SORCE
tsi_path = os.path.join(pf,"sorce_tsi_24hr_l3_full.csv")


fname_AERONET = os.path.join(pf,'aeronet_sites.json')


fname_aerprop_dataset = os.path.join(pf,"TCARS2AERONET_{year}.nc")

with open(fname_AERONET) as txt:
    sites=json.load(txt)                             


# daily integration steps
istep = 0.5 #[h]


years = np.arange(2003,2020)
for year in years:
    # define timeframe
    days = pd.date_range(f"{year}-01-01",f"{year}-12-31")
    if year==2003:
        days = pd.date_range(f"{year}-02-25",f"{year}-12-31")
    for d,day in enumerate(days):
        print(f"  >> {day:%d.%m.%Y}")
        ##########################################################################
        # make ecRad input file

        fn_cams_sfc = os.path.join(ncfile_path,
                                   camsra_fname.format(date=day,levtype='sfc'))
        fn_cams_ml = os.path.join(ncfile_path,
                                  camsra_fname.format(date=day,levtype='ml'))

        sel_lat = np.array([sites[key]['position']['latitude'] for key in sites.keys()])
        sel_lon = np.array([sites[key]['position']['longitude'] for key in sites.keys()])
        # define time, from integration steps
        sel_tim = pd.date_range(day,freq=f'{int(60*istep)}min',periods=int(24/istep))
        # reduce timesteps where sun is definitely not visible.
        # * calculate the solar zenith angle for all time/lat/lon
        # * then we drop all times where sun is not visible
        # * then min/max of this time shows the period of 
        #   maximum sun visibility
        tim,lat = np.meshgrid(sel_tim,sel_lat)
        _,lon = np.meshgrid(sel_tim,sel_lon)
        szen,_ = sp.sun_angles(tim,lat,lon)
        tim_reduced = tim.flatten()[szen.flatten()<100]
        # new time selection based on this selection
        sel_time = pd.date_range(start=np.min(tim_reduced),
                                 end=np.max(tim_reduced),
                                 freq=f'{int(60*istep)}min')  


        # get altitude for each station, repeated for each timestep
        altitudes = np.array([sites[key]['position']['elevation'] for key in sites.keys()])
        sel_alt = np.repeat(altitudes[np.newaxis,:],
                            len(sel_time),axis=0)

        #######################################################################
        # load CAMS RA data and make ecRad input
        # interpolate to dwd station and time selection
        # scale to DWD station surface pressure and altitude
        CamsData = CAMS(fn_cams_sfc,
                        fn_cams_ml,
                        interp_lat=sel_lat,
                        interp_lon=sel_lon,
                        interp_time=sel_time,
                        scale_altitude=sel_alt,
                        scale_sfcpressure=None,
                        combine_latlon=True)


        ap,_ = CamsData.aerosol_optprop([440., 550., 865.],mono=True)


        cams_sfc = CamsData.cams_sfc
        cams_ml = CamsData.cams_ml

        # make ecRad input file
        c2e.make_ecrad_input(date=day,
                             ds_ml=cams_ml,
                             ds_sfc=cams_sfc,
                             outfile_fname='tmp/ecrad_input_aeronet.nc',
                             sorce_tsi_fname=tsi_path)

        # load ecrad input file
        infile = xr.load_dataset('tmp/ecrad_input_aeronet.nc')


        #########################################################################
        # do radiative transfer calculation
        # calculate ecrad only for clear sky times
        eca,ec0 = ecrad.get_ecrad_set(ifile=infile,
                                      reduce_out=False,
                                      ecrad_path=os.path.join(ecrad_path,'bin/'),
                                      opt_prop_file=aerosol_props)

        eca = eca.squeeze()
        ec0 = ec0.squeeze()
        # clean tmp folder
        os.system('rm tmp/*')

        # extract required fluxes
        I_dn_sfc = eca.flux_dn_sw.values[:,-1]
        I0_dn_sfc = ec0.flux_dn_sw.values[:,-1]
        I_up_toa = eca.flux_up_sw.values[:,0]
        I0_up_toa = ec0.flux_up_sw.values[:,0]


        # calculate REari
        REari_sfc = I_dn_sfc - I0_dn_sfc
        REari_toa = I0_up_toa - I_up_toa

        # add station names to dataset 
        names = np.array([key for key in sites.keys()])
        _,names = np.meshgrid(cams_sfc.time.values,names,indexing='ij')


        ap = ap.assign(dict(sites=('column',names.flatten()),
                            REari_sfc=('column',REari_sfc),
                            REari_toa=('column',REari_toa)))

        if d==0:
            ap = ap.assign_coords(column = np.arange(ap.time.shape[0]))
            AP = ap.copy()
        else:
            # append to dataset
            ap = ap.assign_coords(column = np.arange(AP.column[-1]+1,AP.column[-1]+1+ap.time.shape[0]))

            AP = AP.merge(ap,compat='no_conflicts')




        AP.to_netcdf(fname_aerprop_dataset.format(year=year),
                     encoding=dict(column = {'zlib':True},
                                     time = {'zlib':True},
                                     lat = {'zlib':True},
                                     lon = {'zlib':True},
                                     g = {'zlib':True},
                                     ssa = {'zlib':True},
                                     aod = {'zlib':True},
                                     ext = {'zlib':True},
                                  REari_sfc = {'zlib':True},
                                  REari_toa = {'zlib':True}))




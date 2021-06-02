import os
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import configparser
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

# The DWD irradiance observations
fname_DWD_irradiance = os.path.join(pf,'DWD/2015_{station}.nc')


# REari calculated with T-CARS collocated to DWD station
# will be stored in this file
fname_dsTCARS = os.path.join(pf,'REari_TCARS2DWD.nc')


# define timeframe
days = pd.date_range("2015-01-01","2015-12-31")
stations = ld.dwd_stations()
# we drop Zugspitze, since data is not available complete 2015
stations = stations.drop_sel(station='ZG')

# daily integration steps
istep = 0.5 #[h]

def calc_albedo(R,awhite,ablack):
    return (1.-R)*ablack + R*awhite

# variables to store in final dataset
# REari at surface and TOA
variables = ['REari_sfc','REari_toa']


for d,day in enumerate(days):
    print(f"  >> {day:%d.%m.%Y}")
    ##########################################################################
    # make ecRad input file
    
    fn_cams_sfc = os.path.join(ncfile_path,
                               camsra_fname.format(date=day,levtype='sfc'))
    fn_cams_ml = os.path.join(ncfile_path,
                              camsra_fname.format(date=day,levtype='ml'))
                             
    dwdlat = stations.latitude.values
    dwdlon = stations.longitude.values  
    
    # define time, from integration steps
    sel_tim = pd.date_range(day,freq=f'{int(60*istep)}min',periods=int(24/istep))
    # reduce timesteps where sun is definitely not visible.
    # * calculate the solar zenith angle for all time/lat/lon
    # * then we drop all times where sun is not visible
    # * then min/max of this time shows the period of 
    #   maximum sun visibility
    tim,lat = np.meshgrid(sel_tim,dwdlat)
    _,lon = np.meshgrid(sel_tim,dwdlon)
    szen,_ = sp.sun_angles(tim,lat,lon)
    tim_reduced = tim.flatten()[szen.flatten()<100]
    # new time selection based on this selection
    sel_time = pd.date_range(start=np.min(tim_reduced),
                             end=np.max(tim_reduced),
                             freq=f'{int(60*istep)}min')
    
    ## To scale CAMS gridded data to stations, the station altitude and measured pressure
    ## is taken from the DWD dataset
    ## The arrays have to be the shape of the CAMS data after lat/lon interpolation
    
    # get altitude for each station, repeated for each timestep
    dwdalt = np.repeat(stations.height.values[np.newaxis,:],
                       len(sel_time),axis=0)
    # get surface pressure at each station
    for i,station in enumerate(stations.station.values):
        # load DWD station data and select CAMS time
        with xr.open_dataset(fname_DWD_irradiance.format(station=station)) as dwd:
            dwdp0_tmp = dwd.sel(time = sel_time).surf_p_minute.values
            dwdp0_tmp*=1e2 #convert from [hPa] to [Pa]
        if i==0:
            dwdp0 = dwdp0_tmp.copy()
        else:
            dwdp0 = np.vstack((dwdp0,dwdp0_tmp))
    dwdp0 = dwdp0.T
  

    
    #######################################################################
    # load CAMS RA data and make ecRad input
    # interpolate to dwd station and time selection
    # scale to DWD station surface pressure and altitude
    CamsData = CAMS(fn_cams_sfc,
                    fn_cams_ml,
                    interp_lat=dwdlat,
                    interp_lon=dwdlon,
                    interp_time=sel_time,
                    scale_altitude=dwdalt,
                    scale_sfcpressure=dwdp0,
                    combine_latlon=True)


    
    cams_sfc = CamsData.cams_sfc
    cams_ml = CamsData.cams_ml
    
    # make ecRad input file
    c2e.make_ecrad_input(date=day,
                         ds_ml=cams_ml,
                         ds_sfc=cams_sfc,
                         outfile_fname='tmp/ecrad_input_dwd.nc',
                         sorce_tsi_fname=tsi_path)

    # load ecrad input file
    infile = xr.load_dataset('tmp/ecrad_input_dwd.nc')
    
    ######################################################################
    # load lsasaf albedo
    ablack,awhite = ld.load_albedo(day,
                                   coords=(infile.latitude,infile.longitude),
                                   pf=lsasaf_albedo_path)
    
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
    I_dn_direct_sfc = eca.flux_dn_direct_sw.values[:,-1]
    I0_dn_sfc = ec0.flux_dn_sw.values[:,-1]
    I0_dn_direct_sfc = ec0.flux_dn_direct_sw.values[:,-1]
    I_up_toa = eca.flux_up_sw.values[:,0]
    I0_up_toa = ec0.flux_up_sw.values[:,0]
    
    # calculate diff/direct ratio
    with np.errstate(divide='ignore'):
        Raer = (I_dn_sfc - I_dn_direct_sfc) / I_dn_direct_sfc
        Rnoaer = (I0_dn_sfc - I0_dn_direct_sfc) / I0_dn_direct_sfc
    
    # if sun is down, Raer and Rnoaer == nan  (divide by zero)
    # replacing nan with 1 -> only diffuse
    Raer[np.isnan(Raer)] = 1  
    Rnoaer[np.isnan(Rnoaer)] = 1 
    # limit R to [1,0]
    Raer[Raer>1]=1
    Raer[Raer<0]=0
    Rnoaer[Rnoaer>1]=1
    Rnoaer[Rnoaer<0]=0
    
    # calulate albedo
    Aaer = calc_albedo(Raer,awhite,ablack)
    Anoaer = calc_albedo(Rnoaer,awhite,ablack)
    
    # calculate REari
    REari_sfc = (1-Aaer)*I_dn_sfc 
    REari_sfc-= (1-Anoaer)*I0_dn_sfc
    REari_toa = I0_up_toa - I_up_toa
    

    ##########################################################################
    # reshape ecrad output like cams_sfc (time,col)
    newshape = cams_sfc.psfc.shape
    REari_sfc = REari_sfc.reshape(newshape)
    REari_toa = REari_toa.reshape(newshape)

    ###########################################################################
    # daily average
    # integration steps (fraction of day)
    X = (sel_time.hour + sel_time.minute/60)/24.
    
    # do the daily average
    REari_sfc_daily = trapz(REari_sfc,x=X,axis=0)
    REari_toa_daily = trapz(REari_toa,x=X,axis=0)
    
    #####################################################################
    # at first iteration initialyze dataset
    if d==0:
        dummy = np.zeros((len(days),stations.station.size))*np.nan
        dummyvariables={}
        for var in variables:
            dummyvariables.update({var:(('day','station'),dummy.copy())})
        dsTCARS = xr.Dataset(dummyvariables,
                             coords = dict(day=('day',days),
                                           station=('station',stations.station.values)))
    
    
    #####################################################################
    # add to dataset
    dsTCARS.REari_sfc.values[d,:] = REari_sfc_daily
    dsTCARS.REari_toa.values[d,:] = REari_toa_daily

    
encoding = {}
for var in variables:
    encoding.update({var:{'zlib':True}})

dsTCARS.to_netcdf(fname_dsTCARS,
                  encoding=encoding)
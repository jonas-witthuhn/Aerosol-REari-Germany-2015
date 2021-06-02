import os
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import configparser

# custom modules
import modules.load_data as ld
import ecrad_python.camsra2ecrad as c2e
import ecrad_python.run_ecrad as ecrad

# set nice value, since calculation might take a while, and uses cpu alot
os.nice(10)

# paths
config = configparser.ConfigParser()
config.read("ConfigFile.ini")
pf = config['PATHS']['datasets']
ncfile_path = os.path.join(pf,"CAMSRA/")

# path to ecrad installation
# https://github.com/ecmwf/ecrad
ecrad_path = config['PATHS']['ecrad']

# aerosol optical properties of cams aerosol
# https://doi.org/10.24380/jgs8-sc58
# have to be stored in <ecrad_path>/data/
aerosol_props = "aerosol_cams_ifs_optics.nc"

# filename templates
camsra_fname = "cams-ra_{date:%Y-%m-%d}_{levtype}.nc"

# define timeframe
days = pd.date_range("2015-01-01","2015-12-31")

# DWD station metadata
stations = ld.dwd_stations()

print("Run ecRad with and without aerosol:")
newdata=True
dataset_fname = os.path.join(pf,"ecrad_dwd_skill.nc")

for st in stations.station.values:
    print(f"  >> station {st}")
    # select station metadata
    station = stations.sel(station=st)
    lat = station.latitude.values
    lon = station.longitude.values
    
    # clear sky detection mask
    CSD = xr.open_dataset(os.path.join(pf,"CSD/BS_{}.nc".format(st)))
    # DWD irradiance observations
    DWD = xr.open_dataset(os.path.join(pf,"DWD/2015_{}.nc".format(st)))
    
    # loop over every day
    for day in days:
        # select daily data
        csd = CSD.sel(time = day.strftime('%Y-%m-%d'))
        csd = csd.csdc.values # we use the more strict "cloud-free" mask
        dwd = DWD.sel(time = day.strftime('%Y-%m-%d'))

        # skip days with less than 30min of clear sky data
        if np.count_nonzero(csd) < 30:
            continue
        print(f"  >> {day:%d.%m.%Y}")
        ##########################################################################
        # make ecRad input file
        CAMS_ml = xr.open_dataset(os.path.join(ncfile_path,
                                               camsra_fname.format(date=day,
                                                                   levtype='ml')))
        CAMS_sfc = xr.open_dataset(os.path.join(ncfile_path,
                                                camsra_fname.format(date=day,
                                                               levtype='sfc')))

        # interpolate to DWD station
        cams_sfc = CAMS_sfc.interp(lat=lat,lon=lon)
        cams_ml = CAMS_ml.interp(lat=lat,lon=lon)
        
        # select DWD observations at CAMS timesteps
        cams_dwd = dwd.sel(time=cams_sfc.time)
        
        # replace skin temperature and pressure with DWD data
        cams_sfc.psfc.values = cams_dwd.surf_p_minute.values*1e2 #[Pa]
        cams_sfc.tsfc.values = cams_dwd.surf_temp_minute.values+273.15 #[K]
        
        # shift coordinate for albedo at north stations to avoid land - sea confusion
        if stations.sel(tag='north',station=st).site:
            cams_al = CAMS_sfc.interp(lat=lat-1,lon=lon)
            cams_sfc.alnip.values = cams_al.alnip.values
            cams_sfc.alnid.values = cams_al.alnid.values
            cams_sfc.aluvp.values = cams_al.aluvp.values
            cams_sfc.aluvd.values = cams_al.aluvd.values

        # make ecRad input file
        c2e.make_ecrad_input(date=day,
                             ds_ml=cams_ml,
                             ds_sfc=cams_sfc,
                             outfile_fname='tmp/ecrad_input_dwd.nc',
                             sorce_tsi_fname=os.path.join(pf,"sorce_tsi_24hr_l3_full.csv"))
        
        # load ecrad input file
        infile = xr.load_dataset('tmp/ecrad_input_dwd.nc')

        ################################################################################
        # do radiative transfer calculation
        # calculate ecrad only for clear sky times
        times = dwd.time[csd]
        eca,ec0 = ecrad.get_ecrad_set(ifile=infile,
                                      time=times,
                                      reduce_out=True,
                                      ecrad_path=os.path.join(ecrad_path,'bin/'),
                                      opt_prop_file=aerosol_props)
            
        eca = eca.squeeze()
        ec0 = ec0.squeeze()
        mu0 = np.cos(np.deg2rad(eca.sza))
        
        # clean tmp folder
        os.system('rm tmp/ecrad_input_dwd.nc')
        
        ################################################################################
        # store dwd and collocated ecRad simulation
        # glo - global horizontal irradiance
        # dni - direct normal irradiance
        # _obs - DWD observation
        # _aer - ecRad simulation with aerosol
        # _clr - ecRad simulation without aerosol
        
        # offset column dimension to append to complete dataset
        if not newdata:
            N = dsy.col.values[-1]
        else:
            N = 0
        
        ds = xr.Dataset({'glo_obs' : ('col',dwd['global'].values[csd]),
                         'glo_aer' : ('col',eca.flux_dn_sw_sfc.values),
                         'glo_clr' : ('col',ec0.flux_dn_sw_sfc.values),
                         'dni_obs' : ('col',dwd.direct_sc.values[csd]),
                         'dni_aer' : ('col',eca.flux_dn_direct_sw_sfc.values/mu0),
                         'dni_clr' : ('col',ec0.flux_dn_direct_sw_sfc.values/mu0),
                         'mask'    : (('col','tag'),np.ones((len(times),len(station.tag.values))).astype(bool)*station.site.values),
                         'name'    : ('col',[st]*len(times))},
                        coords = {'col' :np.arange(len(times))+N,
                                  'time':('col',times),
                                  'mu0':('col',mu0),
                                  'tag' :station.tag.values})
        # append stats to year
        if newdata:
            dsy = ds.copy()
            newdata = False
        else:
            dsy = xr.combine_by_coords((dsy,ds))

# add compression
encoding={}
for key in dsy.keys():
    encoding.update({key:{'zlib':True}})

# save as netCDF file
dsy.to_netcdf(dataset_fname,encoding=encoding)
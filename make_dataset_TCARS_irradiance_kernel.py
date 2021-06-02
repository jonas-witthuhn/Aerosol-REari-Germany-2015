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

# REari kernel, REari, surface irradiance and aerosol optical props
# will be stored in this file
fname_dsTCARS = os.path.join(pf,'test_TCARS_irradiance_kernel.nc')

# define timeframe
days = pd.date_range("2015-01-01","2015-12-31")

# daily integration steps
istep = 0.5 #[h]

def calc_AE(t1,t2,l1,l2):
    """ calculate Angstrom exponent from aods (t1,t2)
        at two different wavelengths (l1,l2)
    """
    return -np.log(t1/t2)/np.log(l1/l2)

def calc_albedo(DIR,DIF,Adir,Adif):
    """calculate surface albedo using diffuse/direct ratio
    """
    with np.errstate(divide='ignore'):
        S = DIF/DIR
    A = Adir*(1.-S) + Adif*S
    return A

def simulate_reari( infile,
                    perturbation=[],
                    h2o_scale=1,
                    o3_scale=1
                   ):
    """ run ecRad twice, with and without aerosols
        to calculate the REari with the given input file.
        Input parameters can be perturbed / scaled.
        See run_ecrad.get_ecrad_ds for detailed explanation.
    """
    # simulation with aerosols
    OUTF = ecrad.get_ecrad_ds(ifile=infile,
                              outfile = False,
                              reduce_out = True,
                              perturbation=perturbation,
                              aerosol=True,
                              h2o_scale=h2o_scale,
                              o3_scale=o3_scale,
                              co2_scale=1,
                              time=False,
                              ecrad_path = os.path.join(ecrad_path,'bin/'),
                              opfile = aerosol_props,
                              debug=False)
    # spectral net flux
    Fnet_sfc = OUTF.spectral_flux_dn_sw_sfc
    Fnet_sfc-= OUTF.spectral_flux_up_sw_sfc
    Fnet_toa = OUTF.spectral_flux_dn_sw_toa
    Fnet_toa-= OUTF.spectral_flux_up_sw_toa
    
    # simulation without aerosols
    OUTF0 = ecrad.get_ecrad_ds(ifile=infile,
                              outfile = False,
                              reduce_out = True,
                              perturbation=perturbation,
                              aerosol=False,
                              h2o_scale=h2o_scale,
                              o3_scale=o3_scale,
                              co2_scale=1,
                              time=False,
                              ecrad_path = os.path.join(ecrad_path,'bin/'),
                              opfile = aerosol_props,
                              debug=False)
    # spectral net flux
    Fnet0_sfc = OUTF0.spectral_flux_dn_sw_sfc
    Fnet0_sfc-= OUTF0.spectral_flux_up_sw_sfc
    Fnet0_toa = OUTF0.spectral_flux_dn_sw_toa
    Fnet0_toa-= OUTF0.spectral_flux_up_sw_toa
    
    REari_sfc = Fnet_sfc - Fnet0_sfc
    REari_toa = Fnet_toa - Fnet0_toa
    REari = xr.Dataset(dict(TOA = REari_toa,
                            SFC = REari_sfc))
    
    
    mu0 = np.cos(np.deg2rad(OUTF.sza))
    
    # BB downward irradiance components (sfc)
    GLO = OUTF.flux_dn_sw_sfc
    DNI = OUTF.flux_dn_direct_sw_sfc/mu0
    GLO0 = OUTF0.flux_dn_sw_sfc
    DNI0 = OUTF0.flux_dn_direct_sw_sfc/mu0
    I = xr.Dataset(dict(GLO = GLO,
                        GLO0 = GLO0,
                        DNI = DNI,
                        DNI0 = DNI0))
    
    # spectral downward irradiance components (sfc)
    GLO = OUTF.spectral_flux_dn_sw_sfc
    DNI = OUTF.spectral_flux_dn_direct_sw_sfc/mu0
    GLO0 = OUTF0.spectral_flux_dn_sw_sfc
    DNI0 = OUTF0.spectral_flux_dn_direct_sw_sfc/mu0
    Is = xr.Dataset(dict(GLO = GLO,
                        GLO0 = GLO0,
                        DNI = DNI,
                        DNI0 = DNI0))
    return REari,I,Is,mu0


# variables to store in final dataset
# surface irradiance with and without aerosol
variables = ['GLO','DNI','GLO0','DNI0']
# REari at surface and TOA
variables+= ['REari_sfc','REari_toa']
# Aerosol optical properties
variables+= ['AOD','SSA','G','AE']
# Atmospherical parameter
variables+= ['H2O','O3','ALBEDO']
# relative REari kernel [Wm-2] per 1%
variables+= ['kernel_DNI_AOD_sfc',
             'kernel_DNI_SSA_sfc',
             'kernel_DNI_G_sfc',
             'kernel_DNI_AE_sfc',
             'kernel_DNI_H2O_sfc',
             'kernel_DNI_O3_sfc',
             'kernel_DNI_ALBEDO_sfc']
variables+= ['kernel_GLO_AOD_sfc',
             'kernel_GLO_SSA_sfc',
             'kernel_GLO_G_sfc',
             'kernel_GLO_AE_sfc',
             'kernel_GLO_H2O_sfc',
             'kernel_GLO_O3_sfc',
             'kernel_GLO_ALBEDO_sfc']

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
    
    # define time, from integration steps
    sel_tim = pd.date_range(day,freq=f'{int(60*istep)}min',periods=int(24/istep))
    # reduce timesteps where sun is definitely not visible.
    # * calculate the solar zenith angle for all time/lat/lon
    # * then we drop all times where sun is not visible
    # * then min/max of this time shows the period of 
    #   maximum sun visibility
    tim,lat,lon = np.meshgrid(sel_tim,sel_lat,sel_lon)
    szen,_ = sp.sun_angles(tim,lat,lon)
    tim_reduced = tim.flatten()[szen.flatten()<100]
    # new time selection based on this selection
    sel_time = pd.date_range(start=np.min(tim_reduced),
                             end=np.max(tim_reduced),
                             freq=f'{int(60*istep)}min')
    
    

    
    #######################################################################
    # load CAMS RA data and make ecRad input
    # interpolate to reduced grid and time selection
    CamsData = CAMS(fn_cams_sfc,
                    fn_cams_ml,
                    interp_lat=sel_lat,
                    interp_lon=sel_lon,
                    interp_time=sel_time)


    
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
    
    ####################################################################
    # ecRad runs without perturbations
    REari, I, Is ,_ = simulate_reari(infile,
                                perturbation=[],
                                h2o_scale=1,
                                o3_scale=1)
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # calculate bulk aerosol optical properties
    wvls = REari.wvl_sw.values
    AerosolProps,_ = CamsData.aerosol_optprop(wvls,mono=False)
    AP = xr.Dataset(dict(aod = (('time','lat','lon','wvl_sw'),
                                AerosolProps.aod.values.reshape(REari.SFC.shape)),
                         ssa = (('time','lat','lon','wvl_sw'),
                                AerosolProps.ssa.values.reshape(REari.SFC.shape)),
                         g = (('time','lat','lon','wvl_sw'),
                                AerosolProps.g.values.reshape(REari.SFC.shape))
                        ),
                 coords = REari.coords)
    
    # calculate at 440,550,and 865nm to scale BBkernels and calculate AE
    AerosolPropsMono,_ = CamsData.aerosol_optprop([440,550,865],mono=True)
    AP550 = xr.Dataset(dict(aod = (('time','lat','lon'),
                                AerosolPropsMono.aod.values[:,1].reshape(I.GLO.shape)),
                            ssa = (('time','lat','lon'),
                                AerosolPropsMono.ssa.values[:,1].reshape(I.GLO.shape)),
                            g = (('time','lat','lon'),
                                AerosolPropsMono.g.values[:,1].reshape(I.GLO.shape))
                            ),
                        coords = I.coords)                                        
                                                 
    
    # calculate anstrom exponent
    AOD440 = AerosolPropsMono.aod.values[:,0].reshape(I.GLO.shape)
    AOD865 = AerosolPropsMono.aod.values[:,2].reshape(I.GLO.shape)
    AE440865 = calc_AE(AOD440,AOD865,
                       l1=440,l2=865)
    
    
    ####################################################################
    ## kernels scaled to 550 nm
    ####################################################################
    # perturbe AOD by 1%
    perturbation=[('mass_ext_sw_hydrophilic',((slice(None),slice(None),slice(None)),1.)),
                  ('mass_ext_sw_hydrophobic',((slice(None),slice(None)),1.)),
                  ('mass_ext_mono_hydrophilic',((slice(None),slice(None),slice(None)),1.)),
                  ('mass_ext_mono_hydrophobic',((slice(None),slice(None)),1.))
                  ]    
    REari_perturbed,_, I_perturbed,_ = simulate_reari(infile,
                                                    perturbation=perturbation,
                                                    h2o_scale=1,
                                                    o3_scale=1)

    # spectral scaling to 550 nm (Thorsen 2020 Eq.16)
    S = AP.aod/AP550.aod.values[:,:,:,np.newaxis] # add axis dummy for spectral dimension
    
    # calculate kernel scaled to 550nm (Thorsen 2020 Eq.16)
    kernel_AOD = (I_perturbed - Is) * S
    kernel_AOD = kernel_AOD.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # perturbe SSA by 1%
    perturbation = [('ssa_sw_hydrophilic',((slice(None),slice(None),slice(None)),1.)),
                    ('ssa_sw_hydrophobic',((slice(None),slice(None)),1.)),
                    ('ssa_mono_hydrophilic',((slice(None),slice(None),slice(None)),1.)),
                    ('ssa_mono_hydrophobic',((slice(None),slice(None)),1.))
                    ]   
    REari_perturbed,_, I_perturbed,_ = simulate_reari(infile,
                                        perturbation=perturbation,
                                        h2o_scale=1,
                                        o3_scale=1)

    # spectral scaling to 550 nm (Thorsen 2020 Eq.16)
    S = AP.ssa/AP550.ssa.values[:,:,:,np.newaxis] # add axis dummy for spectral dimension
    
    # calculate kernel scaled to 550nm (Thorsen 2020 Eq.16)
    kernel_SSA = (I_perturbed - Is) * S 
    kernel_SSA = kernel_SSA.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # perturbe G by 1%
    perturbation = [('asymmetry_sw_hydrophilic',((slice(None),slice(None),slice(None)),1.)),
                    ('asymmetry_sw_hydrophobic',((slice(None),slice(None)),1.)),
                    ('asymmetry_mono_hydrophilic',((slice(None),slice(None),slice(None)),1.)),
                    ('asymmetry_mono_hydrophobic',((slice(None),slice(None)),1.))
                   ]  
    REari_perturbed,_, I_perturbed,_ = simulate_reari(infile,
                                        perturbation=perturbation,
                                        h2o_scale=1,
                                        o3_scale=1)

    # spectral scaling to 550 nm (Thorsen 2020 Eq.16)
    S = AP.g/AP550.g.values[:,:,:,np.newaxis] # add axis dummy for spectral dimension
    
    # calculate kernel scaled to 550nm (Thorsen 2020 Eq.16)
    kernel_G = (I_perturbed - Is) * S 
    kernel_G = kernel_G.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # perturbe AE by 1%
    perturbation = perturbeAE.perturbe_AE(OPfile=os.path.join(ecrad_path,
                                                              'data/',
                                                              aerosol_props),
                                          pert=1)
    REari_perturbed,_, I_perturbed,_ = simulate_reari(infile,
                                        perturbation=perturbation,
                                        h2o_scale=1,
                                        o3_scale=1)
    
    # not scaled to 550, as AE is considered spectral independend
    kernel_AE = (I_perturbed - Is)
    kernel_AE = kernel_AE.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # perturbe H2O by 1%
    REari_perturbed,_, I_perturbed,_ = simulate_reari(infile,
                                        perturbation=[],
                                        h2o_scale=1.01,
                                        o3_scale=1)
    
    # BB kernel
    kernel_H2O = (I_perturbed - Is)
    kernel_H2O = kernel_H2O.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # perturbe O3 by 1%
    REari_perturbed,_, I_perturbed,_ = simulate_reari(infile,
                                        perturbation=[],
                                        h2o_scale=1,
                                        o3_scale=1.01)
    
    # BB kernel
    kernel_O3 = (I_perturbed - Is)
    kernel_O3 = kernel_O3.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    ####################################################################
    # perturbe surface albedo by 1%
    cams_sfc_perturbed = cams_sfc.copy()
    cams_sfc_perturbed.aluvp.values = cams_sfc_perturbed.aluvp.values*1.01
    cams_sfc_perturbed.aluvd.values = cams_sfc_perturbed.aluvd.values*1.01
    
    c2e.make_ecrad_input(date=day,
                         ds_ml=cams_ml,
                         ds_sfc=cams_sfc_perturbed,
                         outfile_fname='tmp/ecrad_input_dwd.nc',
                         sorce_tsi_fname=tsi_path)

    # load ecrad input file
    infile_perturbed = xr.load_dataset('tmp/ecrad_input_dwd.nc')
    
    REari_perturbed,_, I_perturbed,mu0 = simulate_reari(infile_perturbed,
                                        perturbation=[],
                                        h2o_scale=1,
                                        o3_scale=1)
    
    # BB kernel
    kernel_ALBEDO = (I_perturbed - Is)
    kernel_ALBEDO = kernel_ALBEDO.sum(dim='wvl_sw')
    # clean tmp folder
    os.system('rm tmp/*')
    
    ALBEDO = calc_albedo(DIR=I.DNI.values*mu0,
                         DIF=I.GLO.values-mu0*I.DNI.values,
                         Adir=cams_sfc.aluvp.values,
                         Adif=cams_sfc.aluvd.values)

    
    #####################################################################
    # daily means of atmospheric and aerosol properties
    AE440865 = np.nanmean(AE440865,axis=0)
    AP550 = AP550.mean(dim='time',skipna=True)
    O3 = cams_sfc.o3col.mean(dim='time',skipna=True)
    H2O = cams_sfc.wvcol.mean(dim='time',skipna=True)
    ALBEDO = np.nanmean(ALBEDO,axis=0)
    # daily averages of radiation related parameters 
    # ~daily mean assuming a value of zero if sun is not visible
    REari = REari.sum(dim='wvl_sw').integrate('time',datetime_unit='D')
    I = I.integrate('time',datetime_unit='D')
    kernel_AOD = kernel_AOD.integrate('time',datetime_unit='D')
    kernel_SSA = kernel_SSA.integrate('time',datetime_unit='D')
    kernel_G = kernel_G.integrate('time',datetime_unit='D')
    kernel_AE = kernel_AE.integrate('time',datetime_unit='D')
    kernel_H2O = kernel_H2O.integrate('time',datetime_unit='D')
    kernel_O3 = kernel_O3.integrate('time',datetime_unit='D')
    kernel_ALBEDO = kernel_ALBEDO.integrate('time',datetime_unit='D')
    
    
    #####################################################################
    # at first iteration initialyze dataset
    if d==0:
        dummy = np.zeros((len(days),REari.lat.size,REari.lon.size))*np.nan
        dummyvariables={}
        for var in variables:
            dummyvariables.update({var:(('day','lat','lon'),dummy.copy())})
        dsTCARS = xr.Dataset(dummyvariables,
                             coords = dict(day=('day',days),
                                           lat=('lat',REari.lat),
                                           lon=('lon',REari.lon)))
    
    
    #####################################################################
    # add to dataset
    dsTCARS.REari_sfc.values[d,:,:] = REari.SFC
    dsTCARS.REari_toa.values[d,:,:] = REari.TOA
    dsTCARS.GLO.values[d,:,:] = I.GLO
    dsTCARS.GLO0.values[d,:,:] = I.GLO0
    dsTCARS.DNI.values[d,:,:] = I.DNI
    dsTCARS.DNI0.values[d,:,:] = I.DNI0
    dsTCARS.AOD.values[d,:,:] = AP550.aod
    dsTCARS.SSA.values[d,:,:] = AP550.ssa
    dsTCARS.G.values[d,:,:] = AP550.g
    dsTCARS.AE.values[d,:,:] = AE440865
    dsTCARS.O3.values[d,:,:] = O3
    dsTCARS.H2O.values[d,:,:] = H2O
    dsTCARS.kernel_DNI_AOD_sfc.values[d,:,:] = kernel_AOD.DNI
    dsTCARS.kernel_DNI_SSA_sfc.values[d,:,:] = kernel_SSA.DNI
    dsTCARS.kernel_DNI_G_sfc.values[d,:,:] = kernel_G.DNI
    dsTCARS.kernel_DNI_AE_sfc.values[d,:,:] = kernel_AE.DNI
    dsTCARS.kernel_DNI_H2O_sfc.values[d,:,:] = kernel_H2O.DNI
    dsTCARS.kernel_DNI_O3_sfc.values[d,:,:] = kernel_O3.DNI
    dsTCARS.kernel_DNI_ALBEDO_sfc.values[d,:,:] = kernel_ALBEDO.DNI
    dsTCARS.kernel_GLO_AOD_sfc.values[d,:,:] = kernel_AOD.GLO
    dsTCARS.kernel_GLO_SSA_sfc.values[d,:,:] = kernel_SSA.GLO
    dsTCARS.kernel_GLO_G_sfc.values[d,:,:] = kernel_G.GLO
    dsTCARS.kernel_GLO_AE_sfc.values[d,:,:] = kernel_AE.GLO
    dsTCARS.kernel_GLO_H2O_sfc.values[d,:,:] = kernel_H2O.GLO
    dsTCARS.kernel_GLO_O3_sfc.values[d,:,:] = kernel_O3.GLO
    dsTCARS.kernel_GLO_ALBEDO_sfc.values[d,:,:] = kernel_ALBEDO.GLO
encoding = {}
for var in variables:
    encoding.update({var:{'zlib':True}})

dsTCARS.to_netcdf(fname_dsTCARS,
                  encoding=encoding)
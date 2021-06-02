# Clear sky direct radiative effect of aerosols over Germany

This package includes all python scripts and modules to conduct the analysis done in the paper.
The resulting figures can be found in the *figures* directory.

## get started:
Results as shown in Figures and Tables in the Paper are presented as jupyter notebooks beginning with "results_...".
The required datasets are produced by "make_dataset....py".

To run any of the scripts and notebooks, the raw datasets are required and available from Zenodo: 

So the steps to get started:
1. Clone this repository 
2. Download and extract the datasets folder
3. Edit the paths in "ConfigFile.ini" 
4. Configure a python environment (see below)
5. Try to run the "results_..." notebooks first.
To run the "make_datasets_..." scripts, ecRad is needed:
6. Clone and build the ecrad radiation scheme https://github.com/ecmwf/ecrad 
7. Download the cams aerosol properties file: https://doi.org/10.24380/jgs8-sc58 and store it in the data folder of ecrad

## Dependencies
Use anaconda to create python environment with all required dependencies:
```
 conda env create -f environment.yml
```
This will install python 3.7 with the following modules:
* conda:
    * numpy
    * xarray
    * pandas
    * scipy
    * matplotlib
    * cartopy
    * pyproj
    * h5py
    * netcdf4
    * dask
    * xlsxwriter
    * rpy2
    * jupyter
    * jupyterlab
    * pip
* pip:
    * SkillMetrics (https://pypi.org/project/SkillMetrics/)
    * jstyleson (https://pypi.org/project/jstyleson/)
    * clear_sky_detection (https://github.com/jonas-witthuhn/clear_sky_detection)
    * clear_sky_models (https://github.com/jonas-witthuhn/clear-sky-models)
    * ecrad_python (https://github.com/jonas-witthuhn/ecrad_python)
    * CAMS_aerosol_props (https://github.com/jonas-witthuhn/CAMS_aerosol_props)


## make datasets:
In the following a brief description of the make_datasets scripts are presented.

To calculate all intermediate datasets from scratch, the following order is proposed:
1. [make_dataset_CSD](make_dataset_CSD.py)
2. [make_dataset_lsasaf_albedo](make_dataset_lsasaf_albedo.py)
3. [make_dataset_CSF](make_dataset_CSF.py)
4. [make_dataset_ecrad_skill](make_dataset_ecrad_skill.py)
5. [make_dataset_REari_CSF2DWD](make_dataset_REari_CSF2DWD.py)
6. [make_dataset_AOD_CSM_TCARS](make_dataset_AOD_CSM_TCARS.py)
7. [make_dataset_TCARS_spectral_aerosol_props](make_dataset_TCARS_spectral_aerosol_props.py)
8. [make_dataset_TCARS](make_dataset_TCARS.py)
9. [make_dataset_TCARS2AERONET](make_dataset_TCARS2AERONET.py)
10. [make_dataset_REari_TCARS2DWD](make_dataset_REari_TCARS2DWD.py)

This requires the irradiance observations from DWD (datasets/DWD/), the CAMSRA netcdf data (datasets/CAMSRA/nc/), the LSASAF ALBEDO (datasets/LSASAF) and the AERONET products (datasets/AERONET). The provided datasets include CAMSRA data for 2015 only. 


### make_dataset_AOD_CSM_TCARS.py
Collocated spectral and broadband AOD from TCARS and CSM inversion at DWD stations.
* requires:
    * make_dataset_CSF.py
    * make_dataset_ecrad_skill.py
* produces: TCARS_CSF_AOD_2015.nc
* required for: 
    * results_Fig13.ipynb
    * results_Tab08.ipynb

### make_dataset_CSD.py
Clear sky detection mask (True=clearsky) for all DWD measurements. The clear sky detection is done with the Bright-Sun algorithm for "csds" (free-sun) and "csdc" (cloud-free) conditions.
* requires: -
* produces: CSD/BS_{DWDstation}.nc
* required for: 
    * make_dataset_CSF.py
    * make_dataset_ecrad_skill.py
    * results_Tab02.ipynb

### make_dataset_CSF.py
Conduct clear sky fitting of different clear sky models to the DWD observations, and invert the spectral or broadband AOD.
* requires: 
    * make_dataset_CSD.py
    * make_lsasaf_albedo.py
* produces: CSF/2015_{DWDstation}.nc
* required for: 
    * make_dataset_REari_CSF2DWD.py

### make_dataset_ecrad_skill.py
Run ecrad to simulate clear sky irradiance with CAMSRA input collocated to broadband irradiance measurements from DWD. 
* requires:
    * ecRad
    * make_dataset_CSD.py
* produces: ecrad_dwd_skill.nc
* required for: 
    * results_Tab04_Tab05_Tab06.ipynb

### make_dataset_lsasaf_albedo.py
Extract the surface albedo from the LSASAF surface albedo product, and interpolate to DWD stations.
* requires: - 
* produces: ALBEDO.nc
* required for: 
    * make_dataset_CSF.py
    * make_dataset_REari_CSF2DWD.py

### make_dataset_REari_CSF2DWD.py
Calculate REari from simulated irradiance with and without aerosol with clear sky models from the clear sky fit of the DWD observations.
* requires:
    * make_dataset_CSF.py
    * make_dataset_CSD.py
    * make_dataset_lsasaf_albedo.py
* produces: REari_CSF2DWD.nc
* required for: 
    * results_Fig14_Fig15_Fig16.ipynb
    * results_Tab07_Tab09.ipynb

### make_dataset_REari_TCARS2DWD.py
Simulate REari with TCARS using ecrad with CAMSRA input collocated to DWD stations.
* requires:
    * ecrad
* produces: REari_TCARS2DWD.nc 
* required for: 
    * results_Fig14_Fig15_Fig16.ipynb
    * results_Tab07_Tab09.ipynb

### make_dataset_TCARS_spectral_aerosol_props.py
Calculate spectral aerosol properties (AOD,SSA,G) from CAMSRA mass mixing ratios of the aerosoltypes of the CAMS aerosol model. The calculations are conducted for 2015 over Germany.
* requires: -
* produces:  TCARS_spectral_aerosol_props_2015.nc
* required for: 
    * results_Fig09.ipynb
    * results_FigA1.ipynb

### make_dataset_TCARS.py
Calculate aerosol properties (AOD,SSA,G), direct and global irradiance, REari and REari kernels for 2015 over Germany.
* requires:
    * ecrad
* produces: TCARS.nc
* required for: 
    * results_Fig04_Fig05_FigA2.ipynb
    * results_Fig06_Fig12.ipynb
    * results_Fig08_Fig10_Fig11.ipynb
    * results_Tab07_Tab09.ipynb

### make_dataset_TCARS2AERONET.py
Calculate spectral aerosol properties (AOD,SSA,G) and REari collocated to AERONET stations in Germany for the period from 2003 to 2019. Although, the dataset provided includes CAMSRA data for 2015 only. To calculate for the whole period, one have to download the CAMSRA data from the Copernicus Atmospheric Data Storage. 
* requires:
    * ecrad
* produces: TCARS2AERONET_{year}.nc
* required for: 
    * results_Fig03_Fig07.ipynb

## modules description:
In the following the modules provided in this package are briefly described.
### CAMS_aerosol_props
This module is used to load the CAMS RA data (netcdf). Following features are included:
   * interpolation to any time,lat,lon in the CAMSRA dataset grid
   * scaling the model level data to measured surface pressure and/or altitude
   * calculating AOD,SSA,G at model levels and surface using the CAMS aerosol type mass mixing ratios

### clear_sky_detection
Python version of the Bright-Sun clear sky detection algorithm of measured global+diffuse horizontal irradiance. The original code in matlab can be found at [GitHub](https://github.com/JamieMBright/csd-library)

### clear_sky_models
Python interface to the clear sky models (CSM) library from Jamie Bright coded in R. This is used to calculate the clear sky irradiance at surface. The CSM are used to fill cloud contaminated gaps in the measured irradiance data and estimate REari using a clear sky fitting method. the original code can be fount at [GitHub](https://github.com/JamieMBright/clear-sky-models)

### ecrad_python
Python interface to run the offline version of the [ecRad](https://github.com/ecmwf/ecrad) radiation sheme. Features included:
   * generating ecRad input files from CAMS RA (netcdf) data.
   * interpolating of ecRad input files to required time, lat, lon coordinates
   * reshaping the ecRad output to (time,lat,lon) coordinates
   * perturbations of spectral aerosol properties or atmospheric parameter for sensitivity studies
   
### modules
Several plotting and helperfunctions to load/reshape the data.
   
### retrieve_convert_camsra_camsfc
The scripts in this directory are used to call the Copernicus API to retrieve the CAMSRA data in grib format. The grib files are converted with cdo to netcdf files, which are required for the "CAMS_aerosol_props" module.

### trosat
Functions to calculate sun positon at any time,latitude,longitude




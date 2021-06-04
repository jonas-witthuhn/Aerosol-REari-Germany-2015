## Python modules description:
In the following the modules provided in this package are briefly described.

---

### [CAMS_aerosol_props](https://github.com/jonas-witthuhn/CAMS_aerosol_props)
This module is used to load the CAMS RA data (netcdf). Following features are included:
   * interpolation to any time,lat,lon in the CAMSRA dataset grid
   * scaling the model level data to measured surface pressure and/or altitude
   * calculating AOD,SSA,G at model levels and surface using the CAMS aerosol type mass mixing ratios

---

### [clear_sky_detection](https://github.com/jonas-witthuhn/clear_sky_detection)
Python version of the Bright-Sun clear sky detection algorithm of measured global+diffuse horizontal irradiance. The original code in matlab can be found at [GitHub](https://github.com/JamieMBright/csd-library)

---

### [clear_sky_models](https://github.com/jonas-witthuhn/clear-sky-models)
Python interface to the clear sky models (CSM) library from Jamie Bright coded in R. This is used to calculate the clear sky irradiance at surface. The CSM are used to fill cloud contaminated gaps in the measured irradiance data and estimate REari using a clear sky fitting method. the original code can be found at [GitHub](https://github.com/JamieMBright/clear-sky-models)

---

### [ecrad_python](https://github.com/jonas-witthuhn/ecrad_python)
Python interface to run the offline version of the [ecRad](https://github.com/ecmwf/ecrad) radiation sheme. Features included:
   * generating ecRad input files from CAMS RA (netcdf) data.
   * interpolating of ecRad input files to required time, lat, lon coordinates
   * reshaping the ecRad output to (time,lat,lon) coordinates
   * perturbations of spectral aerosol properties or atmospheric parameter for sensitivity studies
   
---   
   
### modules
Several plotting and helperfunctions to load/reshape the data.

---

### retrieve_convert_camsra_camsfc
The scripts in this directory are used to call the Copernicus API to retrieve the CAMSRA data in grib format. The grib files are converted with cdo to netcdf files, which are required for the "CAMS_aerosol_props" module.

---

### [trosat](https://github.com/hdeneke/trosat-base)
Among others, this package includes functions to calculate sun positon at any time,latitude,longitude. This Analysis exclusively use the trosat.sunpos module for sun position calculation.


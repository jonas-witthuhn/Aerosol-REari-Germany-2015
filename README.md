# Clear sky direct radiative effect of aerosols over Germany

This package includes all python scripts and modules to conduct the analysis done in the paper.
The resulting figures can be found in the *figures* directory.

---

## get started:
Results as shown in Figures and Tables in the Paper are presented as jupyter notebooks beginning with "results_...".
The required intermediate datasets are produced by "make_dataset....py".

To run any of the scripts and notebooks, the raw datasets are required and available from Zenodo: 

So the steps to get started:
1. Clone this repository 
2. Download and extract the datasets folder
3. Edit the paths in [ConfigFile.ini](ConfigFile.ini) 
4. Configure a python environment (see below)
5. Try to run the "results_..." notebooks first.
To run the "make_datasets_..." scripts, ecRad is needed:
6. Clone and build the ecrad radiation scheme https://github.com/ecmwf/ecrad 
7. Download the cams aerosol properties file: https://doi.org/10.24380/jgs8-sc58 and store it in the data folder of ecrad

## Setup Python environment
Use [anaconda](https://www.anaconda.com/) to create python environment with all required dependencies:
```
 conda env create -f environment.yml
 conda acitvate witthuhn2021
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
    * [SkillMetrics](https://pypi.org/project/SkillMetrics/)
    * [jstyleson](https://pypi.org/project/jstyleson/)
    * [clear_sky_detection](https://github.com/jonas-witthuhn/clear_sky_detection)
    * [clear_sky_models](https://github.com/jonas-witthuhn/clear-sky-models)
    * [ecrad_python](https://github.com/jonas-witthuhn/ecrad_python)
    * [CAMS_aerosol_props](https://github.com/jonas-witthuhn/CAMS_aerosol_props)
    * [trosat](https://github.com/hdeneke/trosat-base)

Some modules are created by the authors. For a brief description of the usage of these modules see [PythonModules.md](PythonModules.md).

---

## From scratch:

To follow the analysis of the paper from scratch, the following basis datasets are required:

* **AERONET** - direct sun and inversion products from AERONET 
* **CAMSRA** - CAMS reanalysis model level data for atmospheric and aerosol parameter
* **DWD** - irradiance observations from the German Weather Service (DWD)
* **LSASAF** - LSA SAF surface albedo data 

The data is processed to produce intermediate datasets used for the figure and table creation. Starting from scratch, the following run order is proposed:
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

For a more detailed description of each script see [MakeDataset.md](MakeDataset.md).

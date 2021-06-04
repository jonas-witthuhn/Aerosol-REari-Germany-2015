# make_dataset scripts
In the following a brief description of the make_dataset scripts are presented.

---

### make_dataset_AOD_CSM_TCARS.py
Collocated spectral and broadband AOD from TCARS and CSM inversion at DWD stations.
* requires:
    * (make_dataset_CSF.py)
    * (make_dataset_ecrad_skill.py)
* produces: TCARS_CSF_AOD_2015.nc
* required for: 
    * results_Fig13.ipynb
    * results_Tab08.ipynb

---

### make_dataset_CSD.py
Clear sky detection mask (True=clearsky) for all DWD measurements. The clear sky detection is done with the Bright-Sun algorithm for "csds" (free-sun) and "csdc" (cloud-free) conditions.
* requires: -
* produces: CSD/BS_{DWDstation}.nc
* required for: 
    * make_dataset_CSF.py
    * make_dataset_ecrad_skill.py
    * results_Tab02.ipynb

---

### make_dataset_CSF.py
Conduct clear sky fitting of different clear sky models to the DWD observations, and invert the spectral or broadband AOD.
* requires: 
    * make_dataset_CSD.py
    * make_lsasaf_albedo.py
* produces: CSF/2015_{DWDstation}.nc
* required for: 
    * make_dataset_REari_CSF2DWD.py

---

### make_dataset_ecrad_skill.py
Run ecrad to simulate clear sky irradiance with CAMSRA input collocated to broadband irradiance measurements from DWD. 
* requires:
    * ecRad
    * make_dataset_CSD.py
* produces: ecrad_dwd_skill.nc
* required for: 
    * results_Tab04_Tab05_Tab06.ipynb

---

### make_dataset_lsasaf_albedo.py
Extract the surface albedo from the LSASAF surface albedo product, and interpolate to DWD stations.
* requires: - 
* produces: ALBEDO.nc
* required for: 
    * make_dataset_CSF.py
    * make_dataset_REari_CSF2DWD.py

---

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

---

### make_dataset_REari_TCARS2DWD.py
Simulate REari with TCARS using ecrad with CAMSRA input collocated to DWD stations.
* requires:
    * ecrad
* produces: REari_TCARS2DWD.nc 
* required for: 
    * results_Fig14_Fig15_Fig16.ipynb
    * results_Tab07_Tab09.ipynb

---

### make_dataset_TCARS_spectral_aerosol_props.py
Calculate spectral aerosol properties (AOD,SSA,G) from CAMSRA mass mixing ratios of the aerosoltypes of the CAMS aerosol model. The calculations are conducted for 2015 over Germany.
* requires: -
* produces:  TCARS_spectral_aerosol_props_2015.nc
* required for: 
    * results_Fig09.ipynb
    * results_FigA1.ipynb

---

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

---

### make_dataset_TCARS2AERONET.py
Calculate spectral aerosol properties (AOD,SSA,G) and REari collocated to AERONET stations in Germany for the period from 2003 to 2019. Although, the dataset provided includes CAMSRA data for 2015 only. To calculate for the whole period, one have to download the CAMSRA data from the Copernicus Atmospheric Data Storage. 
* requires:
    * ecrad
* produces: TCARS2AERONET_{year}.nc
* required for: 
    * results_Fig03_Fig07.ipynb

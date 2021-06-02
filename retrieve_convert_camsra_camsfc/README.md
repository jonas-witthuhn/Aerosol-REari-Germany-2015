# retrieve_convert_camsra_camsfc

This is an collection of scripts to retrieve cams-reanalysis and cams-nrealtime with cdsapi and ecmwfapi. Also, there are routines to merge/ convert grib files to netcdf with proper file names.

# Usage:
- see example.py for data download
- see cdo_convert_grib.bash for the merge and conversion to netcdf

# Requirements:
- pandas
- numpy
- cdsapi
- ecmwfapi
- jstyleson

cdsapi and ecmwfapi require an account:
- https://www.ecmwf.int/en/forecasts/access-forecasts/ecmwf-web-api
- https://cds.climate.copernicus.eu/api-how-to


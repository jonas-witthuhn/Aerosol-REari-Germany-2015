import os
import time
import datetime as dt
import concurrent.futures 

from get_cams_grib import get_camsra_date,get_camsfc_date

# Path to store the grib files
pf = "example_data/"
targetfile = os.path.join(pf,'{product}{pfx}_{start:%Y-%m-%d}_{levtype}.grib')
pfx = '_test' # add a tag after the product name (project or something)

# area and grid
area = [56,5.5,46,15.5] # N, W, S, E
grid = [0.25,0.25] # lat/lon grid [deg]

# the downloads will be downloaded day by day
sdates = [dt.date(2019,7,1),
         ]

# initially remove all empty files to clean up failed download attempts
os.system(f"find {pf} -empty -type f -delete")
get_camsra_date(sdates[0],area,grid,targetfile,pfx,['sfc'],True)
# cams-ra
# requires cdsapi
# CDS allowes for 12 requests simultaniously, we issue 13 here so the last one will
# be qeued and imediately issued if one is done, so we dont lose any time here.
with concurrent.futures.ThreadPoolExecutor(max_workers=13) as executor:
    for lev in ['sfc','ml','sfc_met','sfc_rad','sfc_chem']:
        for sd in sdates:
            # for the example, debug is set to True -  no download, just dry run
            executor.submit(get_camsra_date,sd,area,grid,targetfile,pfx,[lev],True)

# cams-fc
# requires ecmwfapi
# Same as for CDS but only 3 requests are allowed simultaniously.
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    for lev in ['sfc','ml','sfc_met','sfc_rad','sfc_chem']:
        for sd in sdates:
            # for the example, debug is set to True -  no download, just dry run
            executor.submit(get_camsfc_date,sd,area,grid,targetfile,pfx,[lev],True)
            time.sleep(5) # issuing the requests to fast will lead them to get blocked

# remove all empty files to clean up failed download attempts
print()
print("Download failed for the following files:")
print(os.system(f"find {pf} -empty -type f"))
print("Try running this script again to fill the gaps.")
print("Cleaning up...")
os.system(f"find {pf} -empty -type f -delete")
print("Done.")

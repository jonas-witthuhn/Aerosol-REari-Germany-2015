#!/bin/bash

date=$1
path="grib/"
tmppath="tmp/"
outpath="nc/"
filetmp=cams-ra_${date}_

### regrid and merge sfc fast and slow request files ####################
echo "regrid sfc data"
# get gridinfo of sfc file
cdo griddes ${path}${filetmp}sfc.grib > "sfcgrid.txt"
# copy sfc (fast) request file to tmp folder
cp ${path}${filetmp}sfc.grib ${tmppath}${filetmp}sfc.grib
fnames=()
# fnames stores all files in the tmp folder
fnames+=(${tmppath}${filetmp}sfc.grib)
# regrid the sfc (slow) request files to the tmp folder
for f in ${path}${filetmp}sfc_*.grib
do
    basename=`basename $f`
    cdo remapbil,"sfcgrid.txt" $f ${tmppath}${basename} # sfcgrid.txt is generated from sfc(fast) request file
    fnames+=(${tmppath}${basename})
done
echo "merge sfc data"
# now all sfc files are on the same grid, so can be merged
tmpsfcfile=${tmppath}${filetmp}sfc_merged.grib
cdo merge ${fnames[@]} ${tmpsfcfile}

### grib to nc ############################################################
echo "convert sfc file to nc"
## use the merged sfc file to convert to nc
cdo -f nc4c -z zip_6 -t cams2cpp.txt copy ${tmpsfcfile} ${outpath}${filetmp}sfc.nc

echo "convert ml file to nc"
cdo -f nc4c -z zip_6 copy ${path}${filetmp}ml.grib ${outpath}${filetmp}ml.nc


### clean up ############################################################
echo "cleanup tmp folder"
rm ${fnames[@]} ${tmpsfcfile}

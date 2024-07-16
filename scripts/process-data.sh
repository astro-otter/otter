#!/bin/bash
set -e

# get the input arguments
INDIR=$1 # path to the unprocessed-data directory in otterdb
OUTDIR=$2 # path to where you want the processed otter jsons stored
FILEDIR=$(dirname -- $0) # the directory where this file is located

echo "Processing data and writing to $OUTDIR..."

# check if the output director exists and if so remove it
if [ -d $OUTDIR ]; then
    echo "$OUTDIR exists already! Removing and refresshing with new data!"
    rm -r $OUTDIR
fi

# make OUTDIR if it doesn't exist
if [ ! -d $OUTDIR ]; then
    echo "Making $OUTDIR because it doesn't exist!"
    mkdir $OUTDIR
fi

###########################################################################
############### FIRST THE UNPROCESSED DATA PEOPLE HAVE SENT ###############
###########################################################################

# first process the base-data which just has some published radio data in it
echo "Processing data from $INDIR/base-data..."
python3 $FILEDIR/base_radio_data_to_otter.py --indir $INDIR/base-data --outdir $OUTDIR

# check that the output files are there
if [[ $(ls $OUTDIR | wc -l) == "1" ]]; then
    echo "Something went wrong, no data found in $OUTDIR!"
    exit 1
fi

# then the data from tde.space
echo "Processing data from $INDIR/tde-1980-2025"
python3 $FILEDIR/tde_dot_space_to_otter.py --indir $INDIR/tde-1980-2025 --outdir $OUTDIR

# then fix the output cause tde.space had some transients in two separate files
python3 $FILEDIR/fix_tde_dot_space_file_bug.py --otterdir $OUTDIR

# then pull data from the "Curated Optical TDE Catalog"
# https://github.com/sjoertvv/manyTDE
wget -P $INDIR https://github.com/sjoertvv/manyTDE/archive/refs/heads/main.zip
unzip $INDIR/main.zip -d $INDIR
rm $INDIR/main.zip
mv $INDIR/manyTDE-main/data/sources $INDIR/curated_optical_tde_catalog
rm -rf $INDIR/manyTDE-main

# and then clean the data from the "Curated Optical TDE Catalog"
python3 $FILEDIR/curated_optical_tde_catalog_to_otter.py --indir $INDIR/curated_optical_tde_catalog --otterdir $OUTDIR

# and then rm the curated_optical_tde_catalog directory since it is already stored
# at the above github
rm -rf $INDIR/curated_optical_tde_catalog

# then the data from Goldtooth et al. (2023)
python3 $FILEDIR/goldtooth_2023_to_otter.py --otterdir $OUTDIR --indir $INDIR/goldtooth_2023

# then the radio data that Noah, Collin, and Kate gathered
python3 $FILEDIR/radio_photometry_to_otter.py --otterdir $OUTDIR --indir $INDIR/radio-data

# then the data for the IR selected TDEs from Masterson et al. (2024)
python3 $FILEDIR/masterson24_to_otter.py --otterdir $OUTDIR --indir $INDIR/masterson24_data

###########################################################################
############### NOW QUERY PUBLIC CATALOGS FOR MORE DATA ###################
###########################################################################
echo "Pulling data from TNS"
python3 $FILEDIR/tns_to_otter.py --otterdir $OUTDIR

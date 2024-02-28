#!/bin/bash
set -e

# get the input arguments
INDIR=$1 # path to the unprocessed-data directory in otterdb
OUTDIR=$2 # path to where you want the processed otter jsons stored
FILEDIR=$(dirname -- $0) # the directory where this file is located

echo "Processing data and writing to $OUTDIR..."

# make OUTDIR if it doesn't exist
if [ ! -d $OUTDIR ]; then
    echo "Making $OUTDIR because it doesn't exist!"
    mkdir $OUTDIR
fi

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

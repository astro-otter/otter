#!/bin/bash

# This script downloads James Guillochan's TDE Catalog and imports
# it to an ArangoDB database

### DOWNLOAD #######################
# change directory and download files
cd ..
git clone git@github.com:astrocatalogs/tde-1980-2025.git

# move into the tde directory
cd tde-1980-2025

### CLEAN #########################
# rename files to get rid of spaces (for ease)
for name in *
do
    mv "$name" "${name// /_}"
done

# add in some brackets for arangodb
for f in *.json
do
    echo -e "[\n$(cat $f)\n]" > $f
done

# remove extra names for ease in queries
for f in *.json
do
    sed -i "$(($(wc -l < $f) - 2))d" $f
    sed -i '3d' $f
done

### IMPORT ########################
# run arango import
for f in *.json
do
    arangoimport --file $f --server.database testing --collection "tdes" --server.username "root@testing" --server.password "password"
done

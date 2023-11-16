#!/bin/bash

### IMPORT ########################
# run arango import
cd ../base
for f in *.json
do
    arangoimport --file $f --server.database "otter" --collection "tdes" --server.username "admin@otter" --server.password "insecure" --define "key=name['default_name']"
done

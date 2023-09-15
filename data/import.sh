#!/bin/bash

### IMPORT ########################
# run arango import
for f in *.json
do
    arangoimport --file $f --server.database "tide" --collection "tdes" --server.username "admin@tide" --server.password "insecure"
done

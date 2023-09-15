#!/bin/bash

# create the database
arangosh --server.endpoint tcp://127.0.0.1:8529 --javascript.execute builddb.js

# import the data
cd ../data
./import.sh

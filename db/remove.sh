#!/bin/bash

arangosh --server.endpoint tcp://127.0.0.1:8529 --javascript.execute cleandb.js

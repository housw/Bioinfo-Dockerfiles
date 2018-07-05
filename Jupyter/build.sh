#!/bin/bash

mkdir -p packages
conda env export --name py3k --file packages/packages.yml
docker build -t 'shengwei/datascience-notebook:latest' .
docker push shengwei/datascience-notebook:latest

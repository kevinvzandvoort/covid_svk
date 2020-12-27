#!/bin/bash

RUN=$1
WD="/nfs/general/covid_svk"
FILE="/nfs/shared/covid_svk/single_run.R"

Rscript ${FILE} ${WD} ${RUN}

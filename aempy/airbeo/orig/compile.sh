#! /bin/bash

# # INTEL OneAPI
# . /opt/intel/oneapi/setvars.sh intel64 ilp64
# ifort -O -g Airbeo.f90 -o iairbeo_orig.x

gfortran -O -g Airbeo.f90 -o gairbeo_orig.x

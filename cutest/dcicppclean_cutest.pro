#!/bin/bash
# Definition for the DCICPP package without NOPE preprocessing

export PACKAGE=dcicppclean
export PACKDIR=dcicppclean
export PACK_PRECISION="double"
export PACK_OBJS=""
export PACKLIBS="-ldcicpp -lnope -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lsuitesparseconfig -lmetis -lopenblas -lgfortran -lgfortranbegin -lpthread -lrt"
export SPECS=""

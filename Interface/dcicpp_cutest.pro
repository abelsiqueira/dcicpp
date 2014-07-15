#!/bin/bash
# Definition for the DCICPP package

export PACKAGE=dcicpp
export PACKDIR=dcicpp
export PACK_PRECISION="double"
export PACK_OBJS=""
export PACKLIBS="-ldcicpp ${HOME}/Libraries/nope/libnope.a -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lsuitesparseconfig -lmetis -lopenblas -lgfortran -lgfortranbegin -lpthread -lrt"
export SPECS=""

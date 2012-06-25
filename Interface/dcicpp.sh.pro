#!/bin/bash
# Definitions for the GENERIC package
# N. Gould, D. Orban & Ph. Toint, November 7th, 2000

# Define a short acronym for the package
export PACK=dcicpp

# Subdirectory of ${CUTER}/common/src/pkg where the package lives
export PACKAGE=dcicpp

# Precision for which the package was written
# Valid values are "single", "double", "single double" and "double single"
export PACK_PRECISION="double"

# Define the name of the object files for the package which must lie in
# ${MYCUTER}/(precision)/bin
export PACKOBJS=""

export PACKLIBS="-ldcicpp -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -I${HOME}/MUMPS_4.10.0/lib -I${HOME}/MUMPS_4.10.0/libseq -I${HOME}/MUMPS_4.10.0/PORT/include -I${HOME}/MUMPS_4.10.0/include ${HOME}/MUMPS_4.10.0/lib/libdmumps.a ${HOME}/MUMPS_4.10.0/lib/libmumps_common.a ${HOME}/MUMPS_4.10.0/lib/libpord.a ${HOME}/MUMPS_4.10.0/libseq/libmpiseq.a -lifcore -limf -lmetis -lmkl_rt -lrt -lsuitesparseconfig /lib64/libpthread.so.0"

# Define the name of the package specification file if any. This possibly
# precision-dependent file must either lie in the current directory or in
# ${CUTER}/common/src/pkg/${PACKAGE}/ )
export SPECS=""

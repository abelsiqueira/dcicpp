# Define DCIDIR before calling this
# Read the INSTALL and define the variables in the $HOME/.bashrc

CXX = icpc

#DCILONG = -DDLONG
DCIVERBOSE = -DVERBOSE
#ASSERT = -DNDEBUG
#PROFILE = -debug -prof-gen=srcpos -prof-dir/home/chico/dcicpp/profile 
PROFILE = 
#OPTIMIZE = -O3
OPTIMIZE = 

#CXXFLAGS = -Werror -Wextra -ggdb $(DCIVERBOSE) $(ASSERT) $(DCILONG) $(PROFILE) -m32
CXXFLAGS = -fexceptions -Wall -Werror $(DCIVERBOSE) $(ASSERT) $(DCILONG) $(BMLONG) $(PROFILE) $(OPTIMIZE)

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

DCILIBNAME = libdcicpp.a

MUMPSDIR=${HOME}/MUMPS_4.10.0
I=-I$(MUMPSDIR)/lib -I$(MUMPSDIR)/libseq -I$(MUMPSDIR)/PORT/include -I$(MUMPSDIR)/include
MUMPSLIB=$(MUMPSDIR)/lib/libdmumps.a $(MUMPSDIR)/lib/libmumps_common.a \
        -L$(MUMPSDIR)/PORD/lib -lpord -L$(MUMPSDIR)/libseq -lmpiseq
#DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis $(MUMPSLIB) -lmkl_rt -lrt -lsuitesparseconfig

#DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -ldmumps -lmumps_common -lpord -lmpiseq -lmkl_rt -lrt -lsuitesparseconfig

DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd ${MUMPSDIR}/lib/libdmumps.a ${MUMPSDIR}/lib/libmumps_common.a ${MUMPSDIR}/lib/libpord.a ${MUMPSDIR}/libseq/libmpiseq.a -lifcore -limf -lmetis -lmkl_rt -lrt -lsuitesparseconfig /lib64/libpthread.so.0

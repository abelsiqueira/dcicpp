# Define DCIDIR before calling this
# Read the INSTALL and define the variables in the $HOME/.bashrc

CXX = g++

#DCILONG = -DLONG
DCIVERBOSE = -DVERBOSE
#ASSERT = -DNDEBUG
#PROFILE = -pg
DEBUG = -ggdb

CXXFLAGS = -Wall -Wextra $(DEBUG) $(DCIVERBOSE) $(ASSERT) $(DCILONG) $(PROFILE) 

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

DCILIBNAME = libdcicpp.a

DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd \
      -lcolamd -lccolamd -lcamd -lsuitesparseconfig -lmetis -ldmumps \
      -lmumps_common -lpord -lmpiseq -lgoto2 -lgfortran -lgfortranbegin \
      -lpthread -lrt

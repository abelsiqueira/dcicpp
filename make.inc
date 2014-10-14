# Define DCIDIR before calling this
# Read the INSTALL and define the variables in the $HOME/.bashrc

CXX = g++
PREFIX = /usr

DCIVERBOSE = -DVERBOSE
#PRINTMATLAB = -DPRINT_MATLAB
#ITERMATLAB = -DITER_MATLAB
#PLOTMATLAB = -DPLOT_MATLAB
#ASSERT = -DNDEBUG
#PROFILE = -pg
DEBUG = -ggdb

CXXFLAGS = -Wall -Wextra $(DEBUG) $(DCIVERBOSE) $(ASSERT) $(PROFILE) $(PRINTMATLAB) \
           $(PLOTMATLAB) $(ITERMATLAB)

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

DCILIBNAME = libdcicpp.a

DCILIB = $(DCIDIR)/lib/$(DCILIBNAME) -lnope -lbasematrices -lcholmod -lamd \
          -lcolamd -lccolamd -lcamd -lsuitesparseconfig -lmetis -lopenblas \
          -lgfortran -lgfortranbegin -lpthread -lrt

I = -I$(DCIDIR)/include/

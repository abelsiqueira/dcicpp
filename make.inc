# Define DCIDIR before calling this
# Read the INSTALL and define the variables in the $HOME/.bashrc

CXX = g++

#DCILONG = -DLONG
DCIVERBOSE = -DVERBOSE
ASSERT = -DNDEBUG
PROFILE = -pg
DEBUG = -ggdb

CXXFLAGS = -Wall -Wextra $(DEBUG) $(DCIVERBOSE) $(ASSERT) $(DCILONG) $(PROFILE) -lrt

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

DCILIBNAME = libdcicpp.a

MUMPSDIR="/home/abel/Libraries/MUMPS_4.10.0"
#I=-I$(MUMPSDIR)/libseq -I$(MUMPSDIR)/PORT/include -I$(MUMPSDIR)/include
#MUMPSLIB=$(MUMPSDIR)/lib/libdmumps.a $(MUMPSDIR)/lib/libmumps_common.a \
        -L$(MUMPSDIR)/PORD/lib -lpord -L$(MUMPSDIR)/libseq -lmpiseq
#DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis $(MUMPSLIB) -lgoto2 -lgfortran -lgfortranbegin -lpthread
DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -ldmumps -lmumps_common -lpord -lmpiseq -lgoto2 -lgfortran -lgfortranbegin -lpthread

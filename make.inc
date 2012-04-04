# Define DCIDIR before calling this
# Read the INSTALL and define the variables in the $HOME/.bashrc

CXX = g++

DCILONG = -DDLONG
#DCIVERBOSE = -DVERBOSE
#ASSERT = -DNDEBUG
#PROFILE = -pg

CXXFLAGS = -Werror -Wextra -ggdb $(DCIVERBOSE) $(ASSERT) $(DCILONG) $(PROFILE) -m32

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

DCILIBNAME = libdcicpp.a

DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -lgoto2 -lgfortran -lgfortranbegin -lpthread

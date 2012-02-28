# Define DCIDIR before calling this
# Read the INSTALL and define the variables in the $HOME/.bashrc

CXX = g++

#DCILONG = -DDLONG
#DCILONG =
#DCIVERBOSE = -DVERBOSE
DCIVERBOSE =
#ASSERT = -DNDEBUG
ASSERT =

CXXFLAGS = -Werror -Wextra -ggdb $(DCIVERBOSE) $(ASSERT) -m32 

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

DCILIBNAME = libdcicpp.a

DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -lgoto2 -lgfortran -lgfortranbegin

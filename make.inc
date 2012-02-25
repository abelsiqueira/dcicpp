# Define DCIDIR before calling this

CXX = g++
#CXXFLAGS = -ggdb

#DCILONG = -DDLONG
#DCILONG =
#DCIVERBOSE = -DVERBOSE
DCIVERBOSE =
#ASSERT = -DNDEBUG
ASSERT =

#CXXFLAGS = -Wall -Werror -Wextra -ggdb $(DCILONG) $(DCIVERBOSE)
#CXXFLAGS = -Wall -Werror -Wextra -ggdb $(DCIVERBOSE) $(ASSERT) -m32 -DLOCALTEST
CXXFLAGS = -Werror -Wextra -ggdb $(DCIVERBOSE) $(ASSERT) -m32 

RM = rm -f
CP = cp -f
AR = ar rv
MV = mv -f

#############################################
##### Adjust this variables accordingly #####
#############################################

##### METIS #####
METIS_PATH = $(HOME)/Libraries/metis
METIS = $(METIS_PATH)/libmetis.a

##### GotoBLAS2 #####
GOTODIR = $(HOME)/Libraries/GotoBLAS2
GOTO = $(GOTODIR)/libgoto2.a $(GOTODIR)/libgoto2.so -lgfortran -lgfortranbegin

##### CHOLMOD, AMD, CAMD, COLAMD, CCOLAMD, UFconfig #####
CHOLMOD_ROOT = $(HOME)/Libraries
AMD = $(CHOLMOD_ROOT)/AMD
CAMD = $(CHOLMOD_ROOT)/CAMD
COLAMD = $(CHOLMOD_ROOT)/COLAMD
CHOLMOD = $(CHOLMOD_ROOT)/CHOLMOD
CCOLAMD = $(CHOLMOD_ROOT)/CCOLAMD
UFCONFIG = $(CHOLMOD_ROOT)/UFconfig
CHOLMODI = -I$(CHOLMOD)/Include -I$(UFCONFIG)
CHOLMODLIB = $(CHOLMOD)/Lib/libcholmod.a $(AMD)/Lib/libamd.a $(CAMD)/Lib/libcamd.a $(COLAMD)/Lib/libcolamd.a $(CCOLAMD)/Lib/libccolamd.a

##### base_matrices #####
BM = $(HOME)/Libraries/base_matrices
BMI = -I$(BM)/Include
BMLIB = $(BM)/Lib/libbasematrices.a

###########################################
##### STOP editing anything past here #####
###########################################
DCILIBNAME = libdcicpp.a
INCLUDES = -I$(DCIDIR)/Include -I$(DCIDIR)/Src $(BMI) -I$(METIS_PATH)/Lib $(CHOLMODI)
LIBS = $(DCIDIR)/Lib/$(DCILIBNAME) $(BMLIB) $(CHOLMODLIB) $(METIS) $(GOTO)
DCII = $(INCLUDES)
DCILIB = $(LIBS)

#DCII = -I$(DCIDIR)/Include -I$(DCIDIR)/Src
#DCILIB = $(DCIDIR)/Lib/$(DCILIBNAME) -lbasematrices -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -lgoto2 -lgfortran -lgfortranbegin

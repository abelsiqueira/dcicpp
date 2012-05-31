#
# Library
#

default: all

DCIDIR = ..
include ../make.inc

OBJ = dcihorz.o dcivert.o porcelli.o restoration.o dcitrust.o dcisteih.o dcisafeguard.o dcibfgs.o method.o solve.o parameters.o interface.o factor.o triplet.o vector.o environment.o 

all: $(OBJ) 
	$(AR) $(DCILIBNAME) $(OBJ)

$(DCILIBNAME): $(OBJ)
	$(AR) $(DCILIBNAME) $(OBJ)

clean:
	$(RM) $(OBJ)

purge: clean
	$(RM) $(DCILIBNAME)

%.o: ../Src/%.cpp 
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(I)

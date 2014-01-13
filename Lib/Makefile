#
# Library
#

default: all

DCIDIR = ..
include ../make.inc

OBJ = dcitangent.o dcinormal.o dciinnernormal.o dciinnertangent.o method.o \
			solve.o parameters.o interface.o factor.o triplet.o vector.o environment.o 

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

#
# Library
#

default: all

DCIDIR = ..
include ../make.inc

OBJS = $(patsubst ../src/%.cpp,%.o,$(wildcard ../src/*.cpp))

all: $(OBJS) $(DCILIBNAME)

$(DCILIBNAME): $(OBJS)
	$(AR) $(DCILIBNAME) $(OBJS)

clean:
	$(RM) $(OBJS)

purge: clean
	$(RM) $(DCILIBNAME)

%.o: ../src/%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(I)

DCIDIR = ..
include ../make.inc

default: all

all: cutest

cutest: dcicppclean dcicpp
	
dcicppclean: dcicppcleanmain_cutest.o
	mkdir -p $(CUTEST)/src/dcicppclean
	$(CP) dcicppclean_cutest.pro $(CUTEST)/packages/$(MYARCH)/double/dcicppclean
	chmod a+x $(CUTEST)/packages/$(MYARCH)/double/dcicppclean
	$(CP) dcicppcleanmain_cutest.cpp $(CUTEST)/src/dcicppclean/dcicppclean_main.cpp
	$(CP) dcicppcleanmain_cutest.o $(CUTEST)/objects/$(MYARCH)/double/dcicppclean_main.o
	$(CP) ../lib/$(DCILIBNAME) $(CUTEST)/objects/$(MYARCH)/double/
	$(CP) makemaster $(CUTEST)/src/dcicppclean/
	$(CP) rundcicppclean $(CUTEST)/bin/

dcicpp: dcicppmain_cutest.o
	mkdir -p $(CUTEST)/src/dcicpp
	$(CP) dcicpp_cutest.pro $(CUTEST)/packages/$(MYARCH)/double/dcicpp
	chmod a+x $(CUTEST)/packages/$(MYARCH)/double/dcicpp
	$(CP) dcicppmain_cutest.cpp $(CUTEST)/src/dcicpp/dcicpp_main.cpp
	$(CP) dcicppmain_cutest.o $(CUTEST)/objects/$(MYARCH)/double/dcicpp_main.o
	$(CP) ../lib/$(DCILIBNAME) $(CUTEST)/objects/$(MYARCH)/double/
	$(CP) makemaster $(CUTEST)/src/dcicpp/
	$(CP) rundcicpp $(CUTEST)/bin/

%.o: %.cpp
	$(CXX) -c -o $@ $< $(DCII) $(CXXFLAGS) $(I) -I$(CUTER)/common/include \
		-I$(CUTEST)/include

clean:
	rm -f *.o

purge: clean

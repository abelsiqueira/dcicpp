DCIDIR = ..
include ../make.inc

default: all

cuter: dcicppmain_cuter.o
	mv dcicppmain_cuter.o dcicppmain.o
	$(CP) dcicpp.sh.pro $(CUTER)/build/prototypes/dcicpp.sh.pro
	sed -f $(MYCUTER)/double/config/script.sed dcicpp.sh.pro > $(MYCUTER)/bin/dcicpp
	chmod a+x $(MYCUTER)/bin/dcicpp
	$(CP) dcicppmain_cuter.cpp $(MYCUTER)/double/bin/dcicppma.cpp
	$(MV) dcicppmain.o $(MYCUTER)/double/bin/dcicppma.o
	$(CP) runcppcuter $(MYCUTER)/bin/
	$(CP) rungdbcuter $(MYCUTER)/bin/
	$(CP) runvalgrindcuter $(MYCUTER)/bin/
	#$(CP) $(GOTOBLASDIR)/libgoto2.a $(MYCUTER)/double/lib/
	#$(CP) $(GOTOBLASDIR)/libgoto2.so $(MYCUTER)/double/lib/
	#$(CP) -t $(MYCUTER)/double/lib/ $(CHOLMODLIB)
	#$(CP) $(BMLIB)/libbasematrices.a $(MYCUTER)/double/lib/

cutest: dcicppmain_cutest.o
	mkdir -p $(CUTEST)/src/dcicpp
	$(CP) dcicpp_cutest.pro $(CUTEST)/packages/$(MYARCH)/double/dcicpp
	chmod a+x $(CUTEST)/packages/$(MYARCH)/double/dcicpp
	$(CP) dcicppmain_cutest.cpp $(CUTEST)/src/dcicpp/dcicpp_main.cpp
	$(CP) dcicppmain_cutest.o $(CUTEST)/objects/$(MYARCH)/double/dcicpp_main.o
	$(CP) ../Lib/$(DCILIBNAME) $(CUTEST)/objects/$(MYARCH)/double/
	$(CP) makemaster $(CUTEST)/src/dcicpp/
	$(CP) rundcicpp $(CUTEST)/bin/

%.o: %.cpp
	$(CXX) -c -o $@ $< $(DCII) $(CXXFLAGS) $(I) -I$(CUTER)/common/include \
		-I$(CUTEST)/include

clean:
	rm -f *.o

purge: clean

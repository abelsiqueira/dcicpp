DCIDIR = ..
include ../make.inc

default: all

all: dcicppmain.o
	$(CP) dcicpp.sh.pro $(CUTER)/build/prototypes/dcicpp.sh.pro
	sed -f $(MYCUTER)/double/config/script.sed dcicpp.sh.pro > $(MYCUTER)/bin/dcicpp
	chmod a+x $(MYCUTER)/bin/dcicpp
	$(CP) dcicppmain.cpp $(MYCUTER)/double/bin/dcicppma.cpp
	$(MV) dcicppmain.o $(MYCUTER)/double/bin/dcicppma.o
	$(CP) runcppcuter $(MYCUTER)/bin/
	$(CP) rungdbcuter $(MYCUTER)/bin/
	$(CP) runvalgrindcuter $(MYCUTER)/bin/
	#$(CP) $(GOTOBLASDIR)/libgoto2.a $(MYCUTER)/double/lib/
	#$(CP) $(GOTOBLASDIR)/libgoto2.so $(MYCUTER)/double/lib/
	#$(CP) -t $(MYCUTER)/double/lib/ $(CHOLMODLIB)
	#$(CP) $(BMLIB)/libbasematrices.a $(MYCUTER)/double/lib/

%.o: %.cpp
	$(CXX) -c -o $@ $< $(DCII) $(CXXFLAGS) -I$(CUTER)/common/include

clean:
	rm -f dcicppmain.o

DCIDIR = ..
include ../make.inc

default: all

#VALGRIND = valgrind
#VALGRINDFLAGS = --leak-check=yes -v
GDB = gdb

TESTSCPP = $(shell ls *test[0-9][0-9].cpp)
TESTS = $(TESTSCPP:.cpp=.tst)
LOCALTEST = localtest.tst
FIXEDTESTS = fixedtest1.tst fixedtest2.tst
EQTESTS = eqtest1.tst eqtest2.tst eqtest3.tst eqtest4.tst eqtest5.tst eqtest6.tst
UNCTESTS = unctest1.tst unctest2.tst
INFEASTESTS = infeastest1.tst infeastest2.tst infeastest3.tst infeastest3e.tst infeastest4.tst infeastest5.tst infeastest6.tst
ARGLINTESTS = arglinb.tst
PIS = points_in_space1.tst points_in_space2.tst
ALLTESTS = $(TESTS) $(OLDTESTS) $(EQTESTS) $(UNCTESTS) $(LOCALTEST) $(INFEASTESTS) $(ARGLINTESTS) $(PIS) $(FIXEDTESTS)

all: library $(TESTS)
	@for i in $(TESTS); do $(VALGRIND) $(VALGRINDFLAGS) ./$$i; done

arglin: library $(ARGLINTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./arglinb.tst

fixedtests: $(FIXEDTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest2.tst

eqtests: $(EQTESTS)
	@for i in $(EQTESTS); do $(VALGRIND) $(VALGRINDFLAGS) ./$$i; done

unctests: $(UNCTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest2.tst

infeastests: $(INFEASTESTS)
	@for i in $(INFEASTESTS); do $(VALGRIND) $(VALGRINDFLAGS) ./$$i; done

pis: $(PIS)
	$(VALGRIND) $(VALGRINDFLAGS) ./points_in_space1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./points_in_space2.tst

localtest: $(LOCALTEST)
	$(VALGRIND) $(VALGRINDFLAGS) ./localtest.tst

library: ../Lib/$(DCILIBNAME)

../Lib/$(DCILIBNAME):
	( cd ..;  make library )


gdb: test14.tst
	$(GDB) ./$<

%.tst: %.o
	$(CXX) -o $@ $< $(DCILIB) $(CXXFLAGS) $(I)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(I)

clean:
	$(RM) $(ALLTESTS) *.o

purge: clean

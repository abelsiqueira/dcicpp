DCIDIR = ..
include ../make.inc

default: all

#VALGRIND = valgrind
#VALGRINDFLAGS = --leak-check=yes -v
GDB = gdb

TESTSCPP = $(shell ls *test[0-9][0-9].cpp)
FIXEDTESTSCPP = $(shell ls fixedtest*.cpp)
EQTESTSCPP = $(shell ls eqtest*.cpp)
UNCTESTSCPP = $(shell ls unctest*.cpp)
INFEASTESTSCPP = $(shell ls infeastest*.cpp)
ALLTESTSCPP = $(shell ls *test*.cpp)
TESTS = $(TESTSCPP:.cpp=.tst)
FIXEDTESTS = $(FIXEDTESTSCPP:.cpp=.tst)
EQTESTS = $(EQTESTSCPP:.cpp=.tst)
UNCTESTS = $(UNCTESTSCPP:.cpp=.tst)
INFEASTESTS = $(INFEASTESTSCPP:.cpp=.tst)
ALLTESTS = $(ALLTESTSCPP:.cpp=.tst)

all: library $(TESTS)
	@for i in $(TESTS); do $(VALGRIND) $(VALGRINDFLAGS) ./$$i; done

fixedtests: $(FIXEDTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest2.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest3.tst

eqtests: $(EQTESTS)
	@for i in $(EQTESTS); do $(VALGRIND) $(VALGRINDFLAGS) ./$$i; done

unctests: $(UNCTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest2.tst

infeastests: $(INFEASTESTS)
	@for i in $(INFEASTESTS); do $(VALGRIND) $(VALGRINDFLAGS) ./$$i; done

library:
	$(MAKE) -C ../lib all

../Lib/$(DCILIBNAME):
	( cd ..;  make library )


gdb: library test01.tst
	$(GDB) ./$(word 2,$^)

%.tst: %.o
	$(CXX) -o $@ $< $(DCILIB) $(CXXFLAGS) $(I)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(I)

clean:
	$(RM) $(ALLTESTS) *.o

purge: clean

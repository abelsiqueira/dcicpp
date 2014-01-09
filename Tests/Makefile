DCIDIR = ..
include ../make.inc

default: all

#VALGRIND = valgrind
#VALGRINDFLAGS = --leak-check=yes -v
GDB = gdb

LOCALTEST = localtest.tst
INTESTS = test1.tst test2.tst test3.tst test4.tst test5.tst test6.tst test7.tst test8.tst test9.tst test10.tst test11.tst test12.tst test13.tst test14.tst test15.tst test16.tst
FIXEDTESTS = fixedtest1.tst fixedtest2.tst fixedtest3.tst
EQTESTS = eqtest1.tst eqtest2.tst eqtest3.tst eqtest4.tst eqtest5.tst test3eq.tst
UNCTESTS = unctest1.tst unctest2.tst
INFEASTESTS = infeastest1.tst infeastest2.tst infeastest3.tst infeastest3e.tst infeastest4.tst infeastest5.tst infeastest6.tst
ARGLINTESTS = arglinb.tst
PIS = points_in_space1.tst points_in_space2.tst
ALLTESTS = $(INTESTS) $(EQTESTS) $(UNCTESTS) $(LOCALTEST) $(INFEASTESTS) $(ARGLINTESTS) $(PIS) $(FIXEDTESTS)

all: library $(INTESTS) 
	$(VALGRIND) $(VALGRINDFLAGS) ./test1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test2.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test3.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test4.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test5.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test6.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test7.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test8.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test9.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test10.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test11.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test12.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test13.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test14.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test15.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./test16.tst

arglin: library $(ARGLINTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./arglinb.tst

fixedtests: $(FIXEDTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest1.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest2.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./fixedtest3.tst 

eqtests: $(EQTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest1.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest2.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest3.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest4.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest5.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./test3eq.tst 

unctests: $(UNCTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest2.tst

infeastests: $(INFEASTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest2.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest3.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest3e.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest4.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest5.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest6.tst

pis: $(PIS)
	$(VALGRIND) $(VALGRINDFLAGS) ./points_in_space1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./points_in_space2.tst

localtest: $(LOCALTEST)
	$(VALGRIND) $(VALGRINDFLAGS) ./localtest.tst

library:
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

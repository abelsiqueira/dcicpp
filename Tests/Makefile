DCIDIR = ..
include ../make.inc

default: all

#VALGRIND = valgrind
#VALGRINDFLAGS = --leak-check=yes -v
GDB = gdb

LOCALTEST = localtest.tst
INTESTS = test1.tst test2.tst test3.tst test4.tst test5.tst test6.tst test7.tst test8.tst test9.tst test10.tst test11.tst test12.tst test13.tst test14.tst
EQTESTS = eqtest1.tst eqtest2.tst eqtest3.tst eqtest4.tst eqtest5.tst
UNCTESTS = unctest1.tst
INFEASTESTS = infeastest1.tst infeastest2.tst infeastest3.tst infeastest3e.tst infeastest4.tst infeastest5.tst
ARGLINTESTS = arglinb.tst
PIS = points_in_space1.tst
ALLTESTS = $(INTESTS) $(EQTESTS) $(UNCTESTS) $(LOCALTEST) $(INFEASTESTS) $(ARGLINTESTS) $(PIS)

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

arglin: library $(ARGLINTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./arglinb.tst

eqtests: $(EQTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest1.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest2.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest3.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest4.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest5.tst 

unctests: $(UNCTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest1.tst

infeastests: $(INFEASTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest1.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest2.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest3.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest3e.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest4.tst
	$(VALGRIND) $(VALGRINDFLAGS) ./infeastest5.tst

pis: $(PIS)
	$(VALGRIND) $(VALGRINDFLAGS) ./points_in_space1.tst

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

DCIDIR = ..
include ../make.inc

default: all

#VALGRIND = valgrind
#VALGRINDFLAGS = --leak-check=yes
VALGRIND = 
VALGRINDFLAGS =
GDB = gdb

LOCALTEST = localtest.tst
INTESTS = test1.tst test2.tst test3.tst test4.tst test5.tst test6.tst test7.tst test8.tst test9.tst test10.tst test11.tst 
EQTESTS = eqtest1.tst eqtest2.tst eqtest3.tst eqtest4.tst eqtest5.tst
UNCTESTS = unctest1.tst
ALLTESTS = $(INTESTS) $(EQTESTS) $(UNCTESTS) $(LOCALTEST)

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

eqtests: $(EQTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest1.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest2.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest3.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest4.tst 
	$(VALGRIND) $(VALGRINDFLAGS) ./eqtest5.tst 

unctests: $(UNCTESTS)
	$(VALGRIND) $(VALGRINDFLAGS) ./unctest1.tst

localtest: $(LOCALTEST)
	$(VALGRIND) $(VALGRINDFLAGS) ./localtest.tst

library:
	( cd ..;  make library )

gdb: unctest1.tst
	$(GDB) ./unctest1.tst

%.tst: %.o
	$(CXX) -o $@ $< $(DCILIB) $(DCII) $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(DCII) $(CXXFLAGS)

clean:
	$(RM) $(ALLTESTS) *.o

purge: clean

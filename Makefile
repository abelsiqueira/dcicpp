default: error

include make.inc

error:
	@echo '------------------------------------------'
	@echo '|          function_handler              |'
	@echo '------------------------------------------'
	@echo '| make all       | make library and test |'
	@echo '| make library   | make library          |'
	@echo '| make test      | test                  |'
	@echo '| make clean     | clean                 |'
	@echo '| make purge     | purge                 |'
	@echo '------------------------------------------'

all: library test

library:
	(cd Lib; make )

gdb: library
	(cd Tests; make gdb)

localtest: library
	(cd Tests; make localtest)

cuter: library
	$(MV) Lib/$(DCILIBNAME) $(MYCUTER)/double/lib/$(DCILIBNAME)
	mkdir -p $(CUTER)/common/src/pkg/dcicpp
	(cd Interface; make cuter)

cutest: library
	#$(MV) Lib/$(DCILIBNAME) $(MYCUTER)/double/lib/$(DCILIBNAME)
	(cd Interface; make cutest)

test: library
	(cd Tests; make )

fixedtests: library
	(cd Tests; make fixedtests)

eqtests: library
	(cd Tests; make eqtests)

unctests: library
	(cd Tests; make unctests)

infeastests: library
	(cd Tests; make infeastests)

arglintests: library
	(cd Tests; make arglin)

clean:
	(cd Lib; make clean )
	(cd Tests; make clean )
	(cd Interface; make clean )

purge:
	(cd Lib; make purge )
	(cd Tests; make purge )
	(cd Interface; make purge )

install: library
	$(CP) Lib/$(DCILIBNAME) $(PREFIX)/lib/
	$(CP) Include/*.h $(PREFIX)/include

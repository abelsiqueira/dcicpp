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
	(cd Interface; make all)

test: library
	(cd Tests; make )

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

purge:
	(cd Lib; make purge )
	(cd Tests; make purge )

default: all

include make.inc

all: library test

library:
	$(MAKE) -C lib all

gdb: library
	$(MAKE) -C tests gdb

cutest: library
	$(MAKE) -C cutest all

test: library
	$(MAKE) -C tests all

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C tests clean
	$(MAKE) -C cutest clean

purge:
	$(MAKE) -C lib purge
	$(MAKE) -C tests purge
	$(MAKE) -C cutest purge

install: library
	$(CP) lib/$(DCILIBNAME) $(PREFIX)/lib/
	$(CP) Include/*.h $(PREFIX)/include

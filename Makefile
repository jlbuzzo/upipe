#!/usr/bin/make
################################################################################


SHELL		:= bash
ECHO		:= echo -e



################################################################################

help:
	$(ECHO) "Usage:\n"
	$(ECHO) "\tmake [options] target [ARGS=\"arg1 arg2 ...\"]"

simple_test:
	$(ECHO) "This is a simple test for arguments:\"$(ARGS)\"."


all:
	$(MAKE) -f scripts/main.mk -C $(PWD) $(ARGS)

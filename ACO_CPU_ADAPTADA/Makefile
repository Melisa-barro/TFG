# Makefile for ACOTSP
OPTIM_FLAGS=-O -lm
WARN_FLAGS=-Wall
CFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS)
CC=gcc
LDLIBS=-lm 

aco: aco.o utilities.o ants.o toymodel.o

all: clean aco

clean:
	@$(RM) *.o aco

aco.o: aco.c

ants.o: ants.c aco.h

utilities.o: utilities.c aco.h

toymodel.o: toymodel.c aco.h


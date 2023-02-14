# Makefile for ACOTSP

# Cuda directories
CUDAPATH = /usr/local/cuda

##########################################################

## C COMPILER OPTIONS ##

# C compiler options:
OPTIM_FLAGS=-O3 -lm
WARN_FLAGS=-Wall
CFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS) -I$(CUDAPATH)/include
CC=gcc

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC=nvcc
NVCCFLAGS = -c -m64 -I$(CUDAPATH)/include

LFLAGS = -m64 -L$(CUDAPATH)/lib -lcuda -lcudart 

all: clean aco

clean:
	@$(RM) *.o aco

aco: ants.o utilities.o toymodel.o aco_main.o
	$(NVCC) $(LFLAGS) aco_main.o ants.o utilities.o toymodel.o -o $@

ants.o: ants.cu aco.h
	$(NVCC) -c $(LFLAGS) $< -o $@

utilities.o: utilities.cu aco.h
	$(NVCC) -c $(LFLAGS) $< -o $@

toymodel.o: toymodel.cu aco.h
	$(NVCC) -c $(LFLAGS) $< -o $@
	
aco_main.o: aco_main.cu aco.h
	$(NVCC) -c $(LFLAGS)  $< -o $@

##########################################################

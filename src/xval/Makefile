# Makefile for hercl/src/xval

ifeq (, $(shell which mpic++))
$(info mpic++ not found,)
ifeq (, $(shell which mpicc))
$(info using gcc instead)
CC = gcc
else
$(info using mpicc instead)
CC = mpicc
endif
else
CC = mpic++
endif

ifeq ($(CC),gcc)
CFLAGS = -g -Wall -O3
else
CFLAGS = -g -DUSE_MPI -Wall -O3 -Wno-write-strings
endif


# Implicit pattern rule
%.o:%.c
	$(CC) $(CFLAGS) -c $<

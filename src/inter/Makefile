# Makefile for hercl/src/inter

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

DEBUG = NO_DEBUG


#paint.o: paint.c Canvas.h Canvas.o
#	$(CC) $(CFLAGS) -c $< Canvas.o

# Implicit pattern rule
%.o:%.c
	$(CC) $(CFLAGS) -D$(DEBUG) -c $<

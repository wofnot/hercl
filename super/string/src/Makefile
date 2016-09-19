# Makefile for hercl/super/data/string/src

CC = gcc
CFLAGS = -g -Wall -O3

string: gen/strcpy gen/strcat gen/strlen gen/idxstr gen/chrstr gen/stridx gen/catstr gen/strchr gen/strrchr gen/strcmp

# Implicit pattern rule
%.o:%.c
	$(CC) $(CFLAGS) -c $<

# implicit pattern rule for gen code
gen/%: %.o mystring.o
	mkdir -p ../gen && $(CC) $(CFLAGS) -lm -o $* $*.o mystring.o && mv $* ../gen

clean:
	rm -f *.o

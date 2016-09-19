# Makefile for hercl/src/base

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
# debug options
#CFLAGS = -gdwarf-3 -Wall -O3
#CFLAGS = -Wall -g -m32
#CFLAGS = -O1 -g -fsanitize=address -fno-omit-frame-pointer
else
CFLAGS = -g -DUSE_MPI -Wall -O3 -Wno-write-strings
endif

# Implicit pattern rule
%.o:%.c
	$(CC) $(CFLAGS) -c $<

# Our target
hercl: hercl.o step.o gradient.o scan_print.o
	$(CC) $(CFLAGS) -o hercl hercl.o gradient.o step.o scan_print.o -lm && mv hercl ../../bin

hercsearch: hercsearch.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o
	$(CC) $(CFLAGS) -lm -o hercsearch hercsearch.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o  && mv hercsearch ../../bin

inter/turtle:../inter/turtle.c intersearch.o interact.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o message.o ../inter/Canvas.o
	cd ../inter && make turtle.o && cd ../base && mv ../inter/turtle.o . && $(CC) $(CFLAGS) -o turtle turtle.o message.o intersearch.o interact.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o  ../inter/Canvas.o -lm && mv turtle ../../bin

# implicit pattern rule for inter code
inter/%:../inter/%.c intersearch.o interact.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o message.o
	cd ../inter && make $*.o && cd ../base && mv ../inter/$*.o . && $(CC) $(CFLAGS) -o $* $*.o message.o intersearch.o interact.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o -lm && mv $* ../../bin

# implicit pattern rule for debug code
debug/%:../inter/%.c debug_interact.o message.o interact.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o
	cd ../inter && make $*.o DEBUG=DEBUG && cd ../base && mv ../inter/$*.o . && $(CC) $(CFLAGS) -o debug_$* $*.o message.o debug_interact.o interact.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o -lm && mv debug_$* ../../bin

# implicit pattern rule for pred code
pred/%:../pred/%.c mystring.o predsearch.o predict.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o
	cd ../pred && make $*.o && cd ../base && mv ../pred/$*.o . && $(CC) $(CFLAGS) -o $* $*.o mystring.o predsearch.o predict.o step.o gradient.o scan_print.o point.o cross_mutate.o eval.o -lm && mv $* ../../bin

# implicit pattern rule for xval code
xval/%:../xval/%.c step.o gradient.o scan_print.o eval.o cross_mutate.o point.o super.o message.o
	cd ../xval && make $*.o && cd ../base && mv ../xval/$*.o . && $(CC) $(CFLAGS) -o $* $*.o step.o gradient.o scan_print.o eval.o cross_mutate.o point.o super.o message.o -lm && mv $* ../../bin

test_message: test_message.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o
	$(CC) $(CFLAGS) -lm -o test_message test_message.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o  && mv test_message ../../bin

master: master.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o
	$(CC) $(CFLAGS) -lm -o master master.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o  && mv master ../../bin

complex: complex.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o
	$(CC) $(CFLAGS) -o complex complex.o scan_print.o step.o gradient.o point.o cross_mutate.o eval.o super.o interact.o inter_hercl.o message.o -lm

society: society.o step.o scan_print.o point.o cross_mutate.o eval.o
	$(CC) $(CFLAGS) -o society society.o step.o scan_print.o point.o cross_mutate.o eval.o -lm

interact.o: interact.c step.h scan_print.h cross_mutate.h eval.h interact.h message.h

intersearch.o: intersearch.c interact.c step.h scan_print.h cross_mutate.h eval.h interact.h

predict.o: predict.c step.h eval.h predict.h

predsearch.o: predsearch.c step.h scan_print.h eval.h predict.h

inter_hercl.o: inter_hercl.c inter_hercl.h step.h eval.h

eval.o: eval.c eval.h step.h scan_print.h cross_mutate.h

gradient.o: gradient.c gradient.h step.h

hercl.o: hercl.c step.h scan_print.h

hercsearch.o: hercsearch.c step.h scan_print.h eval.h super.h interact.h

test_message.o: test_message.c step.h scan_print.h eval.h super.h interact.h message.h

master.o: master.c step.h scan_print.h eval.h super.h interact.h message.h

cross_mutate.o: cross_mutate.c cross_mutate.h step.h scan_print.h point.h

point.o: point.c point.h step.h scan_print.h

scan_print.o: scan_print.c scan_print.h step.h

step.o: step.c step.h gradient.h

super.o: super.c super.h step.h scan_print.h point.h cross_mutate.h eval.h message.h

message.o: message.c message.h eval.h scan_print.h

clean:
	rm -f *.o

clobber:
	cd ../../bin && rm -f hercl hercsearch && cd ../src/base && rm -f *.o

# Makefile for hercl/test

test: quicksort heapsort simqs0 simhs0 simsimqs0 backprop

quicksort: ../bin/hercl
	../bin/hercl -b quicksort.hrc sort.in

heapsort: ../bin/hercl
	../bin/hercl -b heapsort.hrc sort.in

simqs0: ../bin/hercl
	../bin/hercl -b sim.hrc quicksort0.hrc sort.in

simhs0: ../bin/hercl
	../bin/hercl -b sim.hrc heapsort0.hrc sort.in

simsimqs0: ../bin/hercl
	../bin/hercl -b simbig.hrc sim0.hrc quicksort0.hrc sort.in

backprop: ../bin/hercl
	../bin/hercl -b -d 12345 backprop.hrc backprop.cfg ../super/data/raw/xor_.in

../bin/hercl:
	cd ../bin && make hercl && cd ../test

clean:
	cd ../bin && rm -f hercl && cd ../src/base && make clean && cd ../../test

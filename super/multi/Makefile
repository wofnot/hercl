# Makefile for hercl/super/multi

ionosphere: ../../bin/master ../../bin/hercsearch
	mpirun -n 1 ../../bin/master : -n 4 ../../bin/hercsearch -S ionosphere_4.list -j

../../bin/master:
	cd ../../bin && make master && cd ../super/multi

../../bin/hercsearch:
	cd ../../bin && make hercsearch && cd ../super/multi

clean:
	cd ../../bin && rm -f master hercsearch && cd ../src/base && make clean && cd ../../super/multi

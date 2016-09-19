# Makefile for hercl/super/tune

test: xor ionosphere

xor: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/xor.in -u -i ../../lib/hrc/nn16.hrc -a ../../lib/cfg/nn2t1.cfg -t 0

ionosphere: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/ionosphere.in -u -i ../../lib/hrc/nn256.hrc -a ../../lib/cfg/nn4t1.cfg -n

pima: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/pima.in -u -i ../../lib/hrc/nn256.hrc -a ../../lib/cfg/nn4t1.cfg -n -t 0.6

sonar: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/sonar.in -u -i ../../lib/hrc/nn256.hrc -a ../../lib/cfg/nn4t1.cfg -n

../../bin/hercsearch:
	cd ../../bin && make hercsearch && cd ../super/tune

clean:
	cd ../../bin && rm -f hercsearch && cd ../src/base && make clean && cd ../../super/tune


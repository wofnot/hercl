# Makefile for hercl/super/evol

test: xor ionosphere0

xor: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/xor.in -t 0

# train on ionosphere data after splitting into train and test set
ionosphere0: ../../bin/hercsearch
	cd ../data/xval && ./gen_fold.sh ionosphere A && cd ../../evol && ../../bin/hercsearch -s ../data/xval/ionosphere/A/ionosphere0.in -c 4 -n

spirals: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/spirals.in -n -t 0.001

pima: ../../bin/hercsearch
	../../bin/hercsearch -s ../data/raw/pima.in -c 4 -n -t 0.6

../../bin/hercsearch:
	cd ../../bin && make hercsearch && cd ../super/evol

clean:
	cd ../../bin && rm -f hercsearch && cd ../src/base && make clean && cd ../../super/evol

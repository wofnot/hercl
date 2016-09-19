# Makefile for hercl/inter/tune

test: pole_cart pole_cart_hercl

pole_cart_hercl: ../../bin/hercsearch
	../../bin/hercsearch -e ../hrc/pole_cart.hrc -p ../cfg/pole_cart.cfg -u -i ../../lib/hrc/nn16.hrc -a ../../lib/cfg/nn2t1.cfg -n

pole_cart: ../../bin/pole_cart
	../../bin/pole_cart -p ../cfg/pole_cart.cfg -u -i ../../lib/hrc/nn16.hrc -a ../../lib/cfg/nn2t1.cfg -n

majority: ../../bin/majority
	../../bin/majority -u -i ../../lib/hrc/majority.hrc -n 2 -react -t 0.1

../../bin/hercsearch:
	cd ../../bin && make hercsearch && cd ../inter/tune

../../bin/pole_cart:
	cd ../../src/base && make inter/pole_cart && cd ../../inter/tune

clean:
	cd ../../bin && rm -f hercsearch pole_cart && cd ../src/base && make clean && cd ../../inter/tune


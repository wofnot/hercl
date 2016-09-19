# Makefile for hercl/inter/evol

test: pole_cart pole_cart_hercl

pole_cart_hercl: ../../bin/hercsearch
	../../bin/hercsearch -e ../hrc/pole_cart.hrc -p ../cfg/pole_cart.cfg -react -n

pole_cart: ../../bin/pole_cart
	../../bin/pole_cart -p ../cfg/pole_cart.cfg -react -n

double_pole: ../../bin/double_pole
	../../bin/double_pole -p ../cfg/double_pole.cfg -react -n

%:
	cd ../../src/base && make inter/$@ && cd ../../inter/evol

../../bin/hercsearch:
	cd ../../bin && make hercsearch && cd ../inter/evol

../../bin/pole_cart:
	cd ../../src/base && make inter/pole_cart && cd ../../inter/evol

../../bin/double_pole:
	cd ../../src/base && make inter/double_pole && cd ../../inter/evol

../../bin/square_env:
	cd ../../src/base && make inter/square_env && cd ../../inter/evol

../../bin/number_guess:
	cd ../../src/base && make inter/number_guess && cd ../../inter/evol

../../bin/maze_search:
	cd ../../src/base && make inter/maze_search && cd ../../inter/evol

clean:
	cd ../../bin && rm -f hercsearch pole_cart double_pole square_env && cd ../src/base && make clean && cd ../../inter/evol

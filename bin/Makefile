# Makefile for hercl/bin

%:
	if [ -f ../src/base/$*.c ]; then \
            cd ../src/base && make $* ; \
        else \
            if [ -f ../src/cross/$*.c ]; then \
                cd ../src/base && make cross/$* ; \
            else \
                if [ -f ../src/inter/$*.c ]; then \
                    cd ../src/base && make inter/$* ; \
                fi \
            fi \
        fi

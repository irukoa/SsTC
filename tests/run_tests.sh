#!/bin/bash

mv ../Makefile temp_Makefile #Save previous Makefile.

./build_and_compile.sh       #Build and compile.
mpirun -np 1 ./tests.x       #Serial test.
./get_coverage.sh            #Get coverage

#Restore previous Makefile and rebuild.
mv temp_Makefile ../Makefile
(cd ../ && make)

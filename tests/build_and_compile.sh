#!/bin/bash
cp Makefile ../Makefile      #Load testing Makefile and build.
(cd ../ && make)
echo "SsTC built."
#Compile.
mpif90 ./test-drive/src/testdrive.F90 ./suites/suite_l_k_quantities.F90 ./testing_driver.F90 \
                                                       -fprofile-arcs -ftest-coverage \
                                                       -g -fbacktrace --warn-extra --warn-unused --warn-all \
                                                       -Wno-maybe-uninitialized --check=bounds -fstack-check -fimplicit-none \
                                                       -O0 -fopenmp \
                                                       -I../bin/ ../bin/libSsTC.a \
                                                       -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread \
                                                       -o "tests.x"
echo "Test built."

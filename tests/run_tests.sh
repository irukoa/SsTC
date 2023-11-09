#!/bin/bash
#Save previous Makefile.
mv ../Makefile temp_Makefile

echo "Building test suite and compiling..."
echo "-----------------------------------------------"
./build_and_compile.sh > /dev/null
echo "Building done."
echo "-----------------------------------------------"
echo "Running tests (1 MPI process)..."
echo "-----------------------------------------------"
mpirun -np 1 ./tests.x > test_results.dat 2>&1
echo "Running tests (2 MPI processes)..."
echo "-----------------------------------------------"
mpirun -np 2 ./tests.x >> test_results.dat 2>&1
rm toy_model-*
rm s_test.*
rm p_test.*
echo "Done."
echo "-----------------------------------------------"
echo "Getting code coverage..."
echo "-----------------------------------------------"
./get_coverage.sh > coverage.dat 2>&1

echo "Restore previuos setting and rebuild..."
echo "-----------------------------------------------"
mv temp_Makefile ../Makefile
(cd ../ && make > /dev/null)

#Cleanup
rm *.gcda *.gcno
rm *.mod
rm ../*.gcov
shopt -s extglob
cd ../src/obj/
rm !(*.o)
echo "Done."
echo "-----------------------------------------------"

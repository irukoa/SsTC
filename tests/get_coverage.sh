#!/bin/bash
#Retrieve coverage.
(cd ../ && gcov -pb ./src/obj/*.o)
lcov --gcov-tool gcov --capture --directory . --output-file coverage.info
genhtml --output-directory html coverage.info

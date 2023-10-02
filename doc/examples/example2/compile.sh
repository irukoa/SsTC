#!/bin/bash
mpiifort -qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -Ofast example2.F90 -I../../../bin/ ../../../bin/libSsTC.a -o "test.x"

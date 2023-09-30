calculators_jerk.o: $(CALC)/mymod/calculators_jerk.F90 utility.o data_structures.o local_k_quantities.o integrator.o
											 $(F90) $(F90FLAGS) -c $(CALC)/mymod/calculators_jerk.F90 -o "$(OBJ)/calculators_jerk.o"

DEPS += calculators_jerk.o
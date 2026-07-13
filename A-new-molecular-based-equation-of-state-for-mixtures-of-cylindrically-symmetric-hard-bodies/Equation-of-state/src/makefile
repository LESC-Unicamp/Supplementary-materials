.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= cshb.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -c -O3 -ffast-math -funroll-loops

# Object names
OBJS= module_global_variables.o \
	module_substance_properties.o \
	module_phase_equilibria.o \
	module_thermodynamic_processes.o \
	optimization_simplex.o \
	optimization_objective_function.o \
	combining_mixing_rules.o \
	calculation_critical_properties.o \
	calculation_pressure.o \
	calculation_entropy.o \
	calculation_cpotential.o \
	numerical_methods.o \
	numerical_successive_substitution.o \
	numerical_topliss_procedure.o \
	numerical_golden_section_search.o \
	numerical_find_pressure_interval.o \
	numerical_brents_method.o \
	general_routines.o \
	cshb.o \

# Compiling rule
%o: %f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod *.ini $(EXEC_PATH)/*.out

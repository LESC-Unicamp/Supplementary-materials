.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= free-energy-method.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -cpp -fopenmp -c -O3 -ffast-math -funroll-loops

# Object names
OBJS= module_global_variables.o \
	module_quaternion_operations.o \
	module_vector_operations.o \
	module_symmetry_orientfield.o \
	module_overlap_algorithm.o \
	module_overlap_check.o \
	module_initial_configuration.o \
	module_variable_initialization.o \
	module_folders.o \
	module_thermodynamic_path.o \
	subroutine_general.o \
	einsteincrystal.o \

# Compiling rule
%o: %f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -fopenmp -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod $(EXEC_PATH)/*.out

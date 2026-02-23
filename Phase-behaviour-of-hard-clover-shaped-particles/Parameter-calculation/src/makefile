.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= parameters.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -cpp -c -O3 -fopenmp -ffast-math -funroll-loops

# Object names
OBJS= module_global_variables.o \
	module_quaternions.o \
	module_initialization.o \
	module_neighbors.o \
	module_optimization_simplex.o \
	module_vector_operation.o \
	module_image_lists.o \
	module_spherical_harmonics.o \
	module_bond_parameters.o \
	module_structure_factor.o \
	module_order_parameters.o \
	module_radial_distribution_functions.o \
	main.o \

# Compiling rule
%o: %f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -fopenmp -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod $(EXEC_PATH)/*.out

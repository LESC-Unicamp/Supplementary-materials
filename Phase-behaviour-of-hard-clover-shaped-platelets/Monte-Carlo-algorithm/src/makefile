.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= nonconvex.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -cpp -fopenmp -c -O3 -ffast-math -funroll-loops

# Object names
OBJS= module_global_variables.o \
	module_quaternion_operations.o \
	module_vector_operations.o \
	module_latticereduction.o \
	module_folders.o \
	module_linked_lists.o \
	module_overlap_algorithm.o \
	module_overlapcheck_list.o \
	module_overlap_check.o \
	module_initial_configuration.o \
	module_variable_initialization.o \
	subroutine_nematic_order.o \
	subroutine_boxrotation.o \
	subroutine_widom_insertion.o \
	subroutine_general.o \
	nonconvex.o \

# Compiling rule
%o: %f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -fopenmp -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod $(EXEC_PATH)/*.out

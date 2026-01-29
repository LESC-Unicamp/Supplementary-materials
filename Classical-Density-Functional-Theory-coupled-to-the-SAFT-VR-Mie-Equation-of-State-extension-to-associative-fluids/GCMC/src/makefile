.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= gcmc_mie_assoc.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -cpp -c -O3 -ffast-math -funroll-loops

# Object names
OBJS= module_global_variables.o \
	module_quaternion_operations.o \
	module_vector_operations.o \
	subroutine_energy_calculation.o \
	subroutine_force_fields.o \
	subroutine_foursite_model.o \
	subroutine_random_configuration.o \
	subroutine_variable_initialization.o \
	subroutine_widom.o \
	main.o \

# Compiling rule
%.o: %.f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod *.ini $(EXEC_PATH)/*.out

# --------------------------------------------------
# Makefile for quadprog and examples
# --------------------------------------------------

# Files Folder
src = ../src
example = ../examples

# Exec names
exe_F = testf
exe_C = testc

# Compiler
FC = gfortran
CC = gcc
# Flags
W_flags = -Wall -Wextra -Werror
Opt_flags = #-O3 -march=native -flto=8 -fno-fat-lto-objects -fno-strict-aliasing
COM_flags = -c -g -fPIC
F_flags = -ffree-line-length-none #-fcheck=all
C_flags = 
# Dependency
Blas_dir = /opt/apps/blas/3.7.1-gcc-9.1.1
# Valgrind debug
Val_flags = --leak-check=full --show-reachable=yes --track-origins=yes


# Fortran Compiler and Flags
F_Mod = $(FC) $(COM_flags) $(F_flags) $(W_flags) $(Opt_flags)
F_Lin = $(FC) -g $(Opt_flags) -lblas -Wl,-rpath=$(Blas_dir)/lib
# C Compiler and Flags
C_Mod = $(CC) $(COM_flags) $(C_flags) $(W_flags) $(Opt_flags)
C_Lin = $(CC) -g $(Opt_flags) -lgfortran -lblas -Wl,-rpath=$(Blas_dir)/lib

#List of objects to compile
F_list = util.o solve_qp.o quad_prog.o ftest.o
C_list = util.o solve_qp.o quad_prog.o ctest.o

#------------------------
# Build
all: $(F_list) $(C_list)
	$(F_Lin) $(F_list) -o $(exe_F)
	$(C_Lin) $(C_list) -o $(exe_C)
#------------------------

# Generic rules
%.mod: %.o
	@true

%.o: $(src)/%.f
	$(F_Mod) -std=legacy -Wno-extra $<

%.o: $(src)/%.f08
	$(F_Mod) -std=f2008 $<

%.o: $(example)/%.f08
	$(F_Mod) -std=f2008 $<

%.o: $(example)/%.c
	$(C_Mod) $<

# Objects list and dependencies
util.o: $(src)/util.f
solve_qp.o: $(src)/solve_qp.f util.mod
quad_prog.o: $(src)/quad_prog.f08 solve_qp.mod
ftest.o: $(example)/ftest.f08 quad_prog.mod
ctest.o: $(example)/ctest.c util.o solve_qp.o quad_prog.o

#------------------------
# Clean make dir
clean:
	rm *.o *.mod $(exe_F) $(exe_C)

#------------------------
# Run examples
run_C:
	./$(exe_C)
run_F:
	./$(exe_F)

#------------------------
# Run valgrind to debug memory
memdebug_F:
	valgrind $(Val_flags) ./$(exe_F)
memdebug_C:
	valgrind $(Val_flags) ./$(exe_C)

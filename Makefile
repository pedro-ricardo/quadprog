# --------------------------------------------------
# Makefile for ...
# --------------------------------------------------

# Files Folder
src = .

# Exe name
exe_F = testf
exe_C = testc

# Compiler
Fort = gfortran
Cpp = g++
# Flags
W_flags = #-Wall -Wextra -Werror
Opt_flags = #-O3 -march=native -flto=8 -fno-fat-lto-objects -fno-strict-aliasing
COM_flags = -g
For_flags = -ffree-line-length-none #-fcheck=all

# Valgrind debug
Val_flags = --leak-check=full --show-reachable=yes --track-origins=yes

# Concatenate flags and compiler into one command
Fort_Mod = $(Fort) -c $(COM_flags) $(For_flags) $(W_flags) $(Opt_flags)
C_Mod = $(Cpp) -c $(COM_flags) $(W_flags) $(Opt_flags)
Fort_Lin = $(Fort) -g $(Opt_flags)
C_Lin = $(Cpp) -g $(Opt_flags)

#List of objects to compile
For_list = util.o solve_qp.o quad_prog.o ftest.o
C_list = util.o solve_qp.o quad_prog.o control.o

#------------------------
# Build
all: $(For_list) $(C_list)
	$(Fort_Lin) $(For_list) -o $(exe_F) -lblas
	$(C_Lin) $(C_list) -o $(exe_C) -lgfortran -lblas
#------------------------

# Generic rules
%.mod: %.o
	@true

%.o: $(src)/%.f
	$(Fort_Mod) -std=legacy $<

%.o: $(src)/%.f08
	$(Fort_Mod) $<

%.o: $(src)/%.cpp
	$(C_Mod) $<

# Objects list and dependencies
# aind_mod.o: $(src)/aind_mod.f
# solve_qp_compact.o: $(src)/solve_qp_compact.f dpofa_mod.mod dpori_dposl.mod
util.o: $(src)/util.f
solve_qp.o: $(src)/solve_qp.f util.mod
quad_prog.o: $(src)/quad_prog.f08 solve_qp.mod
ftest.o: $(src)/ftest.f08 quad_prog.mod
control.o: $(src)/control.cpp util.o solve_qp.o quad_prog.o

#------------------------
# Clean make dir
clean:
	rm *.o *.mod $(exe_F) $(exe_C)

#------------------------
# Run valgrind to debug memory
mem-debugF:
	valgrind $(Val_flags) ./$(exe_F)
mem-debugC:
	valgrind $(Val_flags) ./$(exe_C)

# --------------------------------------------------
# Makefile for ...
# --------------------------------------------------

# Files Folder
src = .

# Exe name
exe = test

# Compiler
Fort = gfortran
# Flags
W_flags = #-Wall -Wextra -Werror
Opt_flags = #-O3 -march=native -flto=8 -fno-fat-lto-objects -fno-strict-aliasing
All_flags = -g -ffree-line-length-none #-fcheck=all

# Valgrind debug
Val_flags = --leak-check=full --show-reachable=yes --track-origins=yes

# Concatenate flags and compiler into one command
Fort_Mod = $(Fort) -c $(All_flags) $(W_flags) $(Opt_flags)
Fort_Lin = $(Fort) -g $(Opt_flags)

#List of objects to compile
File_list = util.o solve_qp.o quad_prog.o main.o

#------------------------
# Build
all: $(File_list)
	$(Fort_Lin) *.o -o $(exe) -lblas
#------------------------

# Generic rules
%.mod: %.o
	@true

%.o: $(src)/%.f
	$(Fort_Mod) $<

%.o: $(src)/%.f08
	$(Fort_Mod) $< 

# Objects list and dependencies
# aind_mod.o: $(src)/aind_mod.f
# solve_qp_compact.o: $(src)/solve_qp_compact.f dpofa_mod.mod dpori_dposl.mod
util.o: $(src)/util.f
solve_qp.o: $(src)/solve_qp.f util.mod
quad_prog.o: $(src)/quad_prog.f08 solve_qp.mod
main.o: $(src)/main.f08 quad_prog.mod


#------------------------
# Clean make dir
clean:
	rm *.o *.mod $(exe)

#------------------------
# Run valgrind to debug memory
mem-debug:
	valgrind $(Val_flags) ./$(exe)

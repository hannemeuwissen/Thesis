#
# @file Makefile
# @brief Makefile for code to investigate CA-Arnoldi, part of Thesis project in 
# High Performance Computing at Trinity College Dublin.
# @author Hanne Meuwissen (meuwissh@tcd.ie)
# @version 3.0
# @date 2023-05-26
#

CC=gcc
ICC = icc
MPI = mpicc
FLAGS= -g -Wextra -Wall
MKLFLAGS = -mkl
EXECS= caa_gen caa_gen_irreg generate_irr_main distance_subspaces CAA
OBJ = parse.o graph.o sparse.o tsqr_mpi.o bgs.o matrix.o hess.o

all:$(EXECS)

parse.o:parse.c parse.h
	$(CC) -c $< $(FLAGS)

graph.o:graph.c graph.h
	$(CC) -c $< $(FLAGS)

matrix.o:matrix.c matrix.h
	$(ICC) -c $< $(FLAGS)

sparse.o:sparse.c sparse.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS)

bgs.o:bgs.c bgs.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS)

hess.o:hess.c hess.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS)

tsqr_mpi.o:tsqr_mpi.c
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS) -lm

generate_irr_main:generate_irr_main.c graph.o
	$(CC) -o $@ $< graph.o $(FLAGS)

distance_subspaces:distance_subspaces.c
	$(ICC) -o $@ $< $(FLAGS) $(MKLFLAGS) -lm

caa_gen:caa_gen.c $(OBJ)
	$(MPI) -o $@ $< $(OBJ) $(FLAGS) \
        -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -liomp5 -lpthread -lm

caa_gen_irreg:caa_gen_irreg.c $(OBJ)
	$(MPI) -o $@ $< $(OBJ) $(FLAGS) \
        -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -liomp5 -lpthread -lm

CAA:CAA.c $(OBJ)
	$(MPI) -o $@ $< $(OBJ) $(FLAGS) \
        -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -liomp5 -lpthread -lm

.PHONY:clean

clean:
	rm $(EXECS) *.o
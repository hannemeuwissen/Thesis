#
# @file Makefile
# @brief Makefile for thesis at Trinity College Dublin.
# @author Hanne Meuwissen (22307813)
# @version 1.0
# @date 2023-05-26
#

CC=gcc
ICC = icc
MPI = mpicc
FLAGS= -g -Wextra -Wall
MKLFLAGS = -mkl
EXECS= main
OBJ = graph.o sparse.o decomp1d.o tsqr_mpi.o

all:$(EXECS)

decomp1d.o:decomp1d.c
	$(CC) -c $< $(FLAGS)

graph.o:graph.c graph.h
	$(CC) -c $< $(FLAGS)

sparse.o:sparse.c sparse.h
	$(CC) -c $< $(FLAGS)

tsqr_mpi.o:tsqr_mpi.c decomp1d.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS) -lm

main:main.c $(OBJ)
	$(MPI) -o $@ $< $(OBJ) $(FLAGS) $(MKLFLAGS) -fopenmp -lm

.PHONY:clean

clean:
	rm $(EXECS) *.o

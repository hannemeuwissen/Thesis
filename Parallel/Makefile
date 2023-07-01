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
OBJ = parse.o graph.o sparse.o tsqr_mpi.o bgs.o matrix.o hess.o

all:$(EXECS)

parse.o:parse.c parse.h
	$(CC) -c $< $(FLAGS)

matrix.o:matrix.c matrix.h
	$(CC) -c $< $(FLAGS)

graph.o:graph.c graph.h
	$(CC) -c $< $(FLAGS)

sparse.o:sparse.c sparse.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS)

bgs.o:bgs.c bgs.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS)

hess.o:hess.c hess.h
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS)

tsqr_mpi.o:tsqr_mpi.c
	$(ICC) -c $< $(FLAGS) $(MKLFLAGS) -lm

# test:test.c $(OBJ)
# 	$(MPI) -o $@ $< $(OBJ) $(FLAGS) \
#         -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 \
#         -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
#         -liomp5 -lpthread -lm

main:main.c $(OBJ)
	$(MPI) -o $@ $< $(OBJ) $(FLAGS) \
        -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -liomp5 -lpthread -lm

.PHONY:clean

clean:
	rm $(EXECS) *.o

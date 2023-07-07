#
# @file Makefile
# @brief Makefile for code to investigate serial Arnoldi, part of Thesis project in 
# High Performance Computing at Trinity College Dublin.
# @author Hanne Meuwissen (meuwissh@tcd.ie)
# @version 1.0
# @date 2023-05-26
#

CC=gcc
ICC = icc
FLAGS= -g -Wextra -Wall -lm
MKLFLAGS = -mkl
EXECS= main

all:$(EXECS)

graph.o:graph.c graph.h
	$(CC) -c $< $(FLAGS)

sparse.o:sparse.c sparse.h
	$(CC) -c $< $(FLAGS)

main:main.c graph.o sparse.o
	$(ICC) -o $@ $< graph.o sparse.o $(FLAGS) $(MKLFLAGS)

.PHONY:clean

clean:
	rm $(EXECS) *.o

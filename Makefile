#
# @file Makefile
# @brief Makefile for thesis at Trinity College Dublin.
# @author Hanne Meuwissen (22307813)
# @version 1.0
# @date 2023-05-26
#

CC=gcc
ICC = icc
FLAGS= -g -Wextra -Wall
MKLFLAGS = -mkl
EXECS= main

all:$(EXECS)

graph.o:graph.c graph.h
	$(CC) -c $< $(FLAGS)

sparse.o:sparse.c sparse.h

main:main.c graph.o sparse.o
	$(ICC) -o $@ $< graph.o sparse.o $(FLAGS) $(MKLFLAGS)

.PHONY:clean

clean:
	rm $(EXECS) *.o

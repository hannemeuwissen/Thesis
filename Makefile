#
# @file Makefile
# @brief Makefile for thesis at Trinity College Dublin.
# @author Hanne Meuwissen (22307813)
# @version 1.0
# @date 2023-05-26
#

CC=gcc
FLAGS= -g -Wextra -Wall
EXECS= main

all:$(EXECS)

graph.o:graph.c graph.h
	$(CC) -c $< $(FLAGS)

sparse.o:sparse.c sparse.h

main:main.c graph.o sparse.o
	$(CC) -o $@ $< graph.o sparse.o $(FLAGS)

.PHONY:clean

clean:
	rm $(EXECS) *.o

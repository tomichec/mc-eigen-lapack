all: inaEigenLAPACK 
CC=gcc -g
INC=/usr/include/
LIB=-llapacke -llapack -lm

%LAPACK: %LAPACK.c inaEigen.h clancy_rates.h clancy_markov.h
	$(CC) -Wall -I$(INC)/lapacke/ $(LIB) $< -o $@

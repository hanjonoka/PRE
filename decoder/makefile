CXX = gcc
OPT = -O3

GAL_SRC = ../corps_gallois/gallois.c
GAL_HEAD = ../corps_gallois/gallois.h
GRS_H = ../generalized_rs/grs.h
GRS_C = ../generalized_rs/grs.c

decoder: grs galois $(GRS_H) $(GAL_HEAD) multiplicity.c interpolate.c factorize.c main_decoder.c decoder.h
	$(CXX) $(OPT) galois.o grs.o $(GRS_H) $(GAL_HEAD) decoder.h multiplicity.c interpolate.c factorize.c main_decoder.c -o decoder.exe -lm

galois: $(GAL_SRC)
	$(CXX) $(OPT) -c $(GAL_SRC) -o galois.o

grs: galois $(GRS_C) $(GAL_HEAD)
	$(CXX) $(OPT) -c $(GRS_C) -o grs.o
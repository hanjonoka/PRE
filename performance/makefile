CXX = gcc -fopenmp
OPT = -O3

GAL_SRC = ../corps_gallois/gallois.c
GAL_HEAD = ../corps_gallois/gallois.h
GRS_H = ../generalized_rs/grs.h
GRS_C = ../generalized_rs/grs.c
DEC_H = ../decoder/decoder.h
DEC_SRC = ../decoder/multiplicity.c ../decoder/interpolate.c ../decoder/factorize.c ../decoder/decoder.c

DEC_O = multiplicity.o interpolate.o factorize.o decoder.o

test: grs galois decoder $(GRS_H) $(GAL_HEAD) $(DEC_H) test_perf.c
	$(CXX) $(OPT) galois.o grs.o $(DEC_O) $(GRS_H) $(GAL_HEAD) $(DEC_H) test_perf.c -o test_perf.exe -lm

decoder: $(DEC_SRC)
	$(CXX) $(OPT) -c $(DEC_SRC) -lm

galois: $(GAL_SRC)
	$(CXX) $(OPT) -c $(GAL_SRC) -o galois.o

grs: galois $(GRS_C)
	$(CXX) $(OPT) -c $(GRS_C) -o grs.o
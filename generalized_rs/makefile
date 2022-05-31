CXX = gcc
GAL_SRC = ../corps_gallois/gallois.c
GAL_HEAD = ../corps_gallois/gallois.h

rsmain: grs.o gallois.o grsmain.c $(GAL_HEAD) grs.h
	$(CXX) grs.o gallois.o $(GAL_HEAD) grs.h grsmain.c -o grs.exe

grs.o: grs.c
	$(CXX) -c grs.c -o grs.o

gallois.o: $(GAL_SRC)
	$(CXX) -c $(GAL_SRC) -o gallois.o
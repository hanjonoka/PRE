CXX = gcc
GAL_SRC = ../corps_gallois/gallois.cpp
GAL_HEAD = ../corps_gallois/gallois.h

rsmain: grs.o gallois.o grsmain.cpp $(GAL_HEAD) grs.h
	$(CXX) grs.o gallois.o $(GAL_HEAD) grs.h grsmain.cpp -o grs.exe

grs.o: grs.cpp
	$(CXX) -c grs.cpp -o grs.o

gallois.o: $(GAL_SRC)
	$(CXX) -c $(GAL_SRC) -o gallois.o
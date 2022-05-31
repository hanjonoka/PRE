CXX = gcc
GAL_SRC = ../corps_gallois/gallois.c
GAL_HEAD = ../corps_gallois/gallois.h

rsmain: rs.o gallois.o rsmain.c $(GAL_HEAD) rs.h
	$(CXX) rs.o gallois.o $(GAL_HEAD) rs.h rsmain.c -o rs.exe

rs.o: rs.c
	$(CXX) -c rs.c -o rs.o

gallois.o: $(GAL_SRC)
	$(CXX) -c $(GAL_SRC) -o gallois.o
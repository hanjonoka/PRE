CXX = gcc
GAL_SRC = ../corps_gallois/gallois.c
GAL_HEAD = ../corps_gallois/gallois.h

rs: gallois.o rs.c $(GAL_HEAD)
	$(CXX) $(GAL_HEAD)  gallois.o rs.c -o rs.exe

gallois.o: $(GAL_SRC)
	$(CXX) -c $(GAL_SRC) -o gallois.o
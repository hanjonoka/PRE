CXX = gcc
GAL_SRC = ../corps_gallois/gallois.cpp
GAL_HEAD = ../corps_gallois/gallois.h

rsmain: rs.o gallois.o rsmain.cpp $(GAL_HEAD) rs.h
	$(CXX) rs.o gallois.o $(GAL_HEAD) rs.h rsmain.cpp -o rs.exe

rs.o: rs.cpp
	$(CXX) -c rs.cpp -o rs.o

gallois.o: $(GAL_SRC)
	$(CXX) -c $(GAL_SRC) -o gallois.o
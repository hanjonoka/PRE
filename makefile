CXX = gcc
OPT = -O3

GAL_SRC = corps_gallois/gallois.c
GAL_HEAD = corps_gallois/gallois.h
RS_H = reed_solomon/rs.h
RS_C = reed_solomon/rs.c
GRS_H = generalized_rs/grs.h
GRS_C = generalized_rs/grs.c

compute_1: grs rs $(GRS_H) $(RS_H) compute_1.c
	$(CXX) $(OPT) $(GAL_HEAD) $(GRS_H) $(RS_H) rs.o grs.o gallois.o compute_1.c -o compute_1.exe

compare16: grs rs $(GRS_H) $(RS_H) compare16.c
	$(CXX) $(OPT) $(GAL_HEAD) $(GRS_H) $(RS_H) rs.o grs.o gallois.o compare16.c -o compare16.exe -lm

compare_poly: rs $(RS_H) compare_poly.c
	$(CXX) $(OPT) $(GAL_HEAD) $(RS_H) rs.o gallois.o compare_poly.c -o compare_poly.exe

matrix: grs rs $(GRS_H) $(RS_H) matrix.c
	$(CXX) $(OPT) $(GAL_HEAD) $(GRS_H) $(RS_H) rs.o grs.o gallois.o matrix.c -o matrix.exe

compare: grs rs $(GRS_H) $(RS_H) compare.c
	$(CXX) $(OPT) $(GAL_HEAD) $(GRS_H) $(RS_H) rs.o grs.o gallois.o compare.c -o compare.exe -lm

grs: gallois $(GRS_C) $(GAL_HEAD)
	$(CXX) $(OPT) -c $(GRS_C) -o grs.o

rs: gallois $(RS_C) $(GAL_HEAD)
	$(CXX) $(OPT) -c $(RS_C) -o rs.o

gallois: $(GAL_SRC)
	$(CXX) $(OPT) -c $(GAL_SRC) -o gallois.o

clean :
	rm *.o
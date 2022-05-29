CXX = gcc
GAL_SRC = corps_gallois/gallois.c
GAL_HEAD = corps_gallois/gallois.h
RS_H = reed_solomon/rs.h
RS_C = reed_solomon/rs.c
GRS_H = generalized_rs/grs.h
GRS_C = generalized_rs/grs.c

compute_1: grs rs $(RGS_H) $(RS_H) compute_1.c
	$(CXX) $(GAL_HEAD) $(GRS_H) $(RS_H) rs.o grs.o gallois.o compute_1.c -o compute_1.exe

grs: gallois $(GRS_C) $(GAL_HEAD)
	$(CXX) -c $(GRS_C) -o grs.o

rs: gallois $(RS_C) $(GAL_HEAD)
	$(CXX) -c $(RS_C) -o rs.o

gallois: $(GAL_SRC)
	$(CXX) -c $(GAL_SRC) -o gallois.o
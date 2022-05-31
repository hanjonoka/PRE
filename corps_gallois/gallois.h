#ifndef __GALOIS_H__
#define __GALOIS_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct st_galois {
  u_int8_t** elements;
  u_int8_t* add_table;
  u_int8_t* mult_table;
  int n;
  int deg_P;
  u_int8_t* poids;
} galois;


//deprecated
galois* generate_galois_8();
u_int8_t puiss_galois_8(galois* G, u_int8_t x, int n);
u_int8_t pol_of_cart_8(u_int8_t** elts, u_int8_t* e);
void print_galois_8(galois* G);
//addition modulo P = 1 + X^2 + X^3
void mult_galois_8(u_int8_t* a, u_int8_t* b, u_int8_t* res);
//deprecated


galois* generate_galois(u_int8_t* P, int deg_P);
u_int8_t* cart_of_pol(galois* G, int i);
u_int8_t pol_of_cart(galois* G, u_int8_t* e);
void print_galois(galois* G);
u_int8_t puiss_galois(galois* G, u_int8_t x, int n);

#endif
#ifndef __GALLOIS_H__
#define __GALLOIS_H__

#include <stdlib.h>
#include <stdio.h>

typedef struct st_gallois {
  u_int8_t** elements;
  u_int8_t* add_table;
  u_int8_t* mult_table;
} gallois;

gallois* generate_gallois_8();

u_int8_t* cart_of_pol(u_int8_t** elts, int i);

u_int8_t pol_of_cart_8(u_int8_t** elts, u_int8_t* e);

u_int8_t opp_pol_8(u_int8_t i);

void print_gallois_8(gallois* G);

//addition modulo P = 1 + X^2 + X^3
void mult_gallois_8(u_int8_t* a, u_int8_t* b, u_int8_t* res);

#endif
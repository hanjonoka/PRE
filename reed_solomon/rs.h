#ifndef __RS_H__
#define __RS_H__
#include "../corps_gallois/gallois.h"

u_int8_t** encode_rs_8(u_int8_t** msg);
void evaluate_rs_old(u_int8_t** msg, u_int8_t* g, gallois* G, u_int8_t** encoded);


void evaluate_rs(u_int8_t** msg, int msg_len, u_int8_t* g, gallois* G, u_int8_t** encoded);
u_int8_t* get_g(gallois* G, u_int8_t* elts, int deg);
void print_poly(u_int8_t* g, int deg);

#endif
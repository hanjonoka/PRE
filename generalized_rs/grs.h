#ifndef __GRS_H__
#define __GRS_H__
#include "../corps_gallois/gallois.h"

void evaluate_grs_old(u_int8_t** msg, u_int8_t V[7], u_int8_t lambda[7], galois* G, u_int8_t** encoded);
u_int8_t** encode_grs_8(u_int8_t** msg);

void evaluate_grs(u_int8_t** msg, int msg_len, u_int8_t* V, u_int8_t* lambda, galois* G, u_int8_t** encoded);
void evaluate_grs_pol(u_int8_t* msg, int msg_len, u_int8_t* V, u_int8_t* lambda, galois* G, u_int8_t* encoded);
u_int8_t* get_lambda(galois* G, int len, int i);

#endif
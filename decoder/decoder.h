#include "../corps_gallois/gallois.h"

//multiplicity
u_int8_t* generate_multiplicity(galois* G, u_int8_t* received, int k);
int cost(u_int8_t* MM, int height, int width);

//interpolate
int compute_omega(u_int8_t* M, int height, int width, int k);
u_int8_t** interpolate(galois* G, u_int8_t* MM, int height, int width, int K);
void free_poly(u_int8_t** Q, int degx);
u_int8_t** alloc_poly(int degx, int degy);
u_int16_t combi(int n,int k);

//main
void print_reliability_matrix(double* M, int height, int width);
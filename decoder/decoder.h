#include "../corps_gallois/gallois.h"

//multiplicity
u_int8_t* generate_multiplicity(galois* G, u_int8_t* received, int k);
int cost(u_int8_t* MM, int height, int width);
int score(u_int8_t* word, int l, u_int8_t* MM, int height, int width);

//interpolate
int compute_omega(u_int8_t* M, int height, int width, int k);
u_int8_t** interpolate(galois* G, u_int8_t* MM, int height, int width, int K);
void free_poly(u_int8_t** Q, int degx);
u_int8_t** alloc_poly(int degx, int degy);
u_int16_t combi(int n,int k);

//factorize
u_int8_t** Factorize(galois* G, u_int8_t** Q, int degx, int degy, int d, int k, int* n_f);
u_int8_t** Normalize_and_cov(galois* G, u_int8_t** Q, int degx, int degy, u_int8_t f0, int* r_degx, int* r_degy);
u_int8_t** normalize(galois* G, u_int8_t** Q, int degx, int degy, int* r_degx, int* r_degy);

//main
void print_reliability_matrix(double* M, int height, int width);
void print_word(u_int8_t* w, int l);
void print_poly(u_int8_t** Q, int degx, int degy);
void print_poly_x_0(u_int8_t** Q, int degx, int degy);
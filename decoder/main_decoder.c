#include "decoder.h"
#include "../generalized_rs/grs.h"

void print_reliability_matrix(double* M, int height, int width) {
  printf("PI =\n");
  for(int i=0; i<height; i++) {
    for(int j=0; j<width; j++) {
      printf("%lf ", M[i*width + j]);
    }
    printf("\n");
  }
}

void print_multiplicity_matrix(u_int8_t* M, int height, int width) {
  printf("M =\n");
  for(int i=0; i<height; i++) {
    for(int j=0; j<width; j++) {
      printf("%d ", M[i*width + j]);
    }
    printf("\n");
  }
}

void print_poly(u_int8_t** Q, int degx, int degy) {
  printf("Q = ");
  for(int i=0;i<=degx;i++) {
    for(int j=0; j<=degy; j++) {
      if(Q[i][j]!=0) printf(" + %d*x^%d*y^%d", Q[i][j],i,j);
    }
  }
  printf("\n");
}

void print_f(u_int8_t* f, int deg) {
  printf("f(x) = %d",f[0]);
  for(int i=1; i<=deg; i++) {
    printf(" + %d*x^%d", f[i],i);
  }
  printf("\n");
}

void print_list_f(u_int8_t** l_f, int l, int deg) {
  printf("list of f =\n");
  for(int i=0; i<l; i++) {
    printf("\t");
    print_f(l_f[i],deg);
  }
}

void print_poly_x_0(u_int8_t** Q, int degx, int degy) {
  printf("Q = ");
  for(int j=0; j<degy; j++) {
    if(Q[0][j]!=0) printf(" + %d*y^%d", Q[0][j],j);
  }
  printf("\n");
}

void print_word(u_int8_t* w, int l) {
  printf("w = [%d",w[0]);
  for(int i=1; i<l; i++) {
    printf(",%d",w[i]);
  }
  printf("]\n");
}

void print_list_word(u_int8_t** lw, int n, int l) {
  printf("List of words = \n");
  for(int i=0; i<n; i++) {
    printf("\t");
    print_word(lw[i],l);
  }
}

int main() {
  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  u_int8_t* V = get_lambda(G,1);
  u_int8_t* lambda = get_lambda(G,0);
  int k=5;

  u_int8_t msg[5] = {1,2,3,4,5};
  print_f(msg,k-1);

  u_int8_t* received = (u_int8_t*)calloc(G->n-1,sizeof(u_int8_t));
  evaluate_grs_pol(msg,k,V,lambda,G,received);
  print_word(received,7);

  printf("\nMultiplicity Assignement\n");
  u_int8_t* M = generate_multiplicity(G, received, k);
  printf("Result :\n");
  print_multiplicity_matrix(M,G->n-1,G->n);

  printf("\nInterpolation\n");
  u_int8_t** Q = interpolate(G, M, G->n-1, G->n, k);
  int omega = compute_omega(M, G->n-1, G->n, k);
  int L = omega/(k-1);
  int c = cost(M, G->n-1, G->n);
  printf("Result :\n");
  print_poly(Q,c,L);

  printf("\nFactorization\n");
  int n_w;
  int degx_norm, degy_norm;
  u_int8_t** Q_norm = normalize(G, Q, c, L, &degx_norm, &degy_norm);
  print_poly(Q_norm, degx_norm, degy_norm);
  u_int8_t** l_w = Factorize(G, Q_norm, degx_norm, degy_norm, 0, k, &n_w);
  print_list_f(l_w, n_w, k-1);
  return 0;
}
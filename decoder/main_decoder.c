#include "decoder.h"

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
    for(int j=0; j<degy; j++) {
      if(Q[i][j]!=0) printf(" + %d*x^%d*y^%d", Q[i][j],i,j);
    }
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

int main() {
  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  int k=5;

  u_int8_t received[7] = {6,5,4,5,0,4,6};
  print_word(received,7);

  u_int8_t* M = generate_multiplicity(G, received, k);
  print_multiplicity_matrix(M,G->n-1,G->n);

  u_int8_t** Q = interpolate(G, M, G->n-1, G->n, k);
  int omega = compute_omega(M, G->n-1, G->n, k);
  int L = omega/(k-1);
  int c = cost(M, G->n-1, G->n);
  print_poly(Q,c,L);
  return 0;
}
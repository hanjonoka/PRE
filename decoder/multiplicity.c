#include "decoder.h"
#include "math.h"

//matrice de taille n*q, q=G->n, n=q-1
double* generate_reliabitity_matrix(galois* G, u_int8_t* received) {
  double p_error = 0.3;

  double* matrix = calloc(G->n*G->n, sizeof(double));
  for(int i=0; i<G->n-1; i++) {
    for(int j=0; j<G->n; j++) {
      u_int8_t d = G->distances[received[i]*G->n + j];
      matrix[i*G->n + j] = pow(p_error,d) * pow(1.0 - p_error,G->deg_P - d);
    }
  }

  return matrix;
}

int cost(u_int8_t* MM, int height, int width) {
  int c=0;
  for(int i=0; i<width*height; i++) {
    c += (MM[i]*(MM[i]+1))/2;
  }
  return c;
}

int nb_monome(int a, int b, int omega) { //nb de monome de wdeg(a,b) au plus omega
  int i=1;
  int j=0;
  int nb=0;

  while(i*a < omega) {
    while(j*b < omega) {
      nb++;
      j++;
    }
    j=0;
    i++;
  }

  return nb;
}

int compute_omega(u_int8_t* M, int height, int width, int k) {
  int c = cost(M, height, width);
  int omega=0;
  while(nb_monome(1,k-1,omega)<c) {
    omega++;
  }
  return omega;
}

u_int8_t* generate_multiplicity(galois* G, u_int8_t* received, int k) {

  double* PI = generate_reliabitity_matrix(G, received); print_reliability_matrix(PI, G->n-1, G->n);
  u_int8_t* MM = calloc(G->n*G->n, sizeof(u_int8_t));

  while(compute_omega(MM, G->n-1, G->n, k) < 20) { //permet de s'arreter quand le cout devient trop important.

    int i_max = 0;
    int j_max = 0;
    double ratio_max = PI[i_max*G->n + j_max]/(MM[i_max*G->n + j_max]+2);

    for(int i=0; i<G->n; i++) { //calcul de i_max, j_max qui maximisent le rapport gain/cout
      for(int j=0; j<G->n; j++) {
        double r = PI[i*G->n + j]/(MM[i*G->n + j]+2);
        if(r>ratio_max) {
          i_max=i;
          j_max=j;
          ratio_max=r;
        }
      }
    }

    MM[i_max*G->n + j_max]++;
  }

  return MM;
}
#include "decoder.h"
#include "math.h"

#define pi 3.14159265358979323846

//matrice de taille n*q, q=G->n, n=q-1
double* generate_reliabitity_matrix_hard(galois* G, u_int8_t* received) {
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

double* generate_reliabitity_matrix_soft(galois* G, double* received, double sigma) {
  double* matrix = calloc(G->n*G->n, sizeof(double));
  for(int i=0; i<G->n-1; i++) {
    double mu = received[i];
    double s=0;
    for(int j=0; j<G->n; j++) {
      matrix[i*G->n + j] = 1/(sigma*sqrt(2*pi))*exp(-pow(j-mu,2)/(2*pow(sigma,2)));
      s+=matrix[i*G->n + j];
    }
    for(int j=0; j<G->n; j++) {
      matrix[i*G->n + j] = matrix[i*G->n + j]/s;
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

int score(u_int8_t* word, int l, u_int8_t* MM, int height, int width) {
  int score=0;
  for(int i=0; i<l; i++) {
    score+=MM[i*width + word[i]];
  }
  return score;
}

int nb_monome(int a, int b, int omega) { //nb de monome de wdeg(a,b) au plus omega
  if(b==0 ^ a==0) return omega;
  int i=0;
  int j=0;
  int nb=0;

  while(i*a <= omega) {
    while(i*a + j*b <= omega) {
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

u_int8_t* generate_multiplicity_greedy(galois* G, u_int8_t* received, int k) {

  double* PI = generate_reliabitity_matrix_hard(G, received); print_reliability_matrix(PI, G->n-1, G->n);
  u_int8_t* MM = calloc((G->n-1)*G->n, sizeof(u_int8_t));

  while(compute_omega(MM, G->n-1, G->n, k) < 3*G->n) { //permet de s'arreter quand le cout devient trop important.

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

u_int8_t* generate_multiplicity_proportional(galois* G, double* PI, int max_mult, int k) {
  double max_pi=PI[0];
  int i_max_pi=0;
  for(int i=0; i<G->n*G->n-1; i++) {
    if(PI[i]>max_pi){
      i_max_pi=i;
      max_pi=PI[i];
    }
  }
  double lambda = ((double)max_mult)/max_pi;

  u_int8_t* M = calloc((G->n-1)*G->n, sizeof(u_int8_t));
  for(int i=0; i<(G->n-1)*G->n; i++) {
    M[i]=lambda*PI[i];
  }

  int omega = compute_omega(M,G->n-1,G->n,k);
  int i_max = 0;
  int j_max = 0;
  while(compute_omega(M, G->n-1, G->n, k) == omega) { //permet de s'arreter quand omega change.
    double ratio_max = PI[i_max*G->n + j_max]/(M[i_max*G->n + j_max]+2);

    for(int i=0; i<G->n; i++) { //calcul de i_max, j_max qui maximisent le rapport gain/cout
      for(int j=0; j<G->n; j++) {
        double r = PI[i*G->n + j]/(M[i*G->n + j]+2);
        if(r>ratio_max) {
          i_max=i;
          j_max=j;
          ratio_max=r;
        }
      }
    }

    M[i_max*G->n + j_max]++;
  }
  M[i_max*G->n + j_max]--;

  return M;
}
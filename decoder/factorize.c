#include "decoder.h"
#include "../generalized_rs/grs.h"

u_int8_t Q_of_Fq_Fq(galois* G, u_int8_t** Q, int degx, int degy, u_int8_t X, u_int8_t Y) {
  u_int8_t res = 0;

  for(int  i=0; i<degx; i++) {
    for(int j=0; j<degy; j++) {
      u_int8_t tX = puiss_galois(G,X,i);
      u_int8_t tY = puiss_galois(G,Y,j);
      u_int8_t r = G->mult_table[tX*G->n + tY];
      r = G->mult_table[Q[i][j]*G->n + r];
      res = G->add_table[r*G->n + res];
    }
  }
  return res;
}

u_int8_t Q_of_0_Fq(galois* G, u_int8_t** Q, int degx, int degy, u_int8_t Y) {
  u_int8_t res=0;
  for(int j=0; j<degy; j++) {
    u_int8_t tY = puiss_galois(G,Y,j);
    u_int8_t r = G->mult_table[Q[0][j]*G->n + r];
    res = G->add_table[r*G->n + res];
  }
  return res;
}

//returns the list of all roots f0_i in Fq of Q(0,y) by testing all elements in Fq
u_int8_t* find_y_roots(galois* G, u_int8_t** Q, int degx, int degy, int* n_roots) {
  *n_roots = 0;
  u_int8_t* l_roots = calloc(G->n, sizeof(u_int8_t));
  for(int i=0; i<G->n; i++) {
    if(Q_of_0_Fq(G, Q, degx, degy, i) == 0) {
      l_roots[*n_roots]=i;
      (*n_roots)++;
    }
  }
  return l_roots;
}

//Normalizes Q(x,xy+f0), f0 in Fq. returns a new polynome.
u_int8_t** Normalize_and_cov(galois* G, u_int8_t** Q, int degx, int degy, u_int8_t f0, int* r_degx, int* r_degy) {
  //change of variable
  int newdegx = degx + degy;
  int newdegy = degy;

  u_int8_t** cov_Q = alloc_poly(newdegx, newdegy);
  //q_ij*x^i*(xy+f0)^j = \SUM_k=0^j [combi(k,j)*q_ij*f0^(j-k) * x^(i+k) * y^(k)]
  for(int i=0; i<=degx; i++) {
    for(int j=0; j<=degy; j++) {
      //binome de newton
      for(int k=0; k<=j; k++) {
        u_int8_t coef = G->mult_table[Q[i][j]*G->n + puiss_galois(G,f0,j-k)]; //q_ij * f0^(j-k)
        coef = mult_n_galois(G, coef, combi(j, k));
        cov_Q[i+k][k] = G->add_table[cov_Q[i+k][k]*G->n + coef];
      }
    }
  }

  //normalization
  //trouver m max tq x^m|Q(x,y)
  //trouver m max tq pour i<=m et pour tout j Q[i][j]=0
  int flag = 0;
  int m=-1;
  for(int i=0; i<=newdegx && !flag; i++) {
    m++;
    for(int j=0; j<=newdegy; j++) {
      if(cov_Q[i][j]!=0) flag=1;
    }
  }

  //division par x^m
  *r_degx = newdegx-m;
  *r_degy = newdegy;
  u_int8_t** r_Q = alloc_poly(*r_degx, *r_degy);
  for(int i=0; i<=*r_degx; i++) {
    for(int j=0; j<=*r_degy; j++) {
      r_Q[i][j] = cov_Q[i+m][j];
    }
  }

  free_poly(cov_Q,newdegx);
  return r_Q;
}

//Finds all y-roots f0 in Fq[x] of Q of degree less than k
u_int8_t** Factorize(galois* G, u_int8_t** Q, int degx, int degy, int i, int k, int* n_f) {
  if(i>=k) {
    *n_f = 0;
    return calloc(*n_f, sizeof(u_int8_t*));
  }

  int n_roots;
  u_int8_t* l_roots = find_y_roots(G, Q, degx, degy, &n_roots);
  u_int8_t** l_f = NULL;
  *n_f = 0;
  for(int i=0; i<n_roots; i++) {

    //calcul liste polynomes commençant par l_roots[i]
    int n_p;
    u_int8_t** l_p = Factorize(G, , , , i+1, k, &n_p);
    *n_f = *n_f + n_p;
    l_f = realloc(*n_f, sizeof(u_int8_t*));

    //concaténation de l_roots[i] et p pour tout p dans l_p
  }
}
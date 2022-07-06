#include "decoder.h"
#include "../generalized_rs/grs.h"


u_int8_t*** init_polys(int degx, int degy) {

  u_int8_t*** l_polys = (u_int8_t***) calloc(degy+1,sizeof(u_int8_t**));
  for(int i=0; i<degy+1; i++) {
    l_polys[i] = alloc_poly(degx, degy);
    l_polys[i][0][i] = 1;
  }

  return l_polys;
}

u_int8_t** alloc_poly(int degx, int degy) {
  u_int8_t** Q = calloc(degx+1, sizeof(u_int8_t*));
  for(int i=0; i<=degx; i++) {
    Q[i] = calloc(degy+1, sizeof(u_int8_t));
  }
  return Q;
}

void free_poly(u_int8_t** Q, int degx) {
  for(int i=0; i<=degx; i++) {
    free(Q[i]);
  }
  free(Q);
}

int weighted_degree(u_int8_t** Q, int degx, int degy, int k) {
  int wd_max = 0;
  for(int i=0; i<=degx; i++) {
    for(int j=0; j<=degy; j++) {
      if(Q[i][j]!=0) {
        int wd=1*i + (k-1)*j;
        if(wd>wd_max) wd_max=wd;
      }
    }
  }
  return wd_max;
}

int deg_y(u_int8_t** Q, int degx, int degy) {
  int deg_max = 0;
  for(int i=0; i<degx; i++) {
    for(int j=0; j<degy; j++) {
      if(Q[i][j]!=0) {
        if(j>deg_max) deg_max=j;
      }
    }
  }
  return deg_max;
}

int greater_than(u_int8_t** Q1, u_int8_t** Q2, int degx, int degy, int k) {
  int wd1 = weighted_degree(Q1, degx, degy, k);
  int wd2 = weighted_degree(Q2, degx, degy, k);

  if(wd1==wd2) {
    return deg_y(Q1, degx, degy)>deg_y(Q2, degx, degy);
  } else {
    return wd1>wd2;
  }
}

u_int16_t combi(int n,int k)
{
    u_int16_t ans=1;
    k=k>n-k?n-k:k;
    int j=1;
    for(;j<=k;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return ans;
}

u_int8_t D(galois* G, u_int8_t** Q, int degx, int degy, int u, int v, u_int8_t alpha, u_int8_t beta) {
  u_int8_t res = 0;
  for(int i=u; i<=degx; i++) {
    for(int j=v; j<=degy; j++) {
      u_int8_t a = G->mult_table[puiss_galois(G,alpha,i-u)*G->n + puiss_galois(G,beta,j-v)];
      a = G->mult_table[Q[i][j]*G->n + a];
      a = mult_n_galois(G, a, combi(i,u));
      a = mult_n_galois(G, a, combi(j,v));
      res = G->add_table[res*G->n + a];
    }
  }
  return res;
}

u_int8_t** interpolate(galois* G, u_int8_t* MM, int height, int width, int K) {
  //init
  int omega = compute_omega(MM, height, width, K);
  int L = omega/(K-1);
  int c = cost(MM, height, width);
  printf("omega=%d, L=%d, c=%d, k=%d\n",omega,L,c,K);

  u_int8_t* discrepancy = calloc(L,sizeof(u_int8_t));
  u_int8_t*** l_polys = init_polys(c, L);
  // for(int k=0; k<=L; k++) {
  //   print_poly(l_polys[k],c,L);
  // }

  u_int8_t* V = get_lambda(G,1);


  //parcours de la matrice de multiplicité
  for(int i=0; i<G->n-1; i++) {
    for(int j=0; j<G->n; j++) {
      u_int8_t alpha = V[i];
      u_int8_t beta = j;


      //generation des contraintes en alpha,beta
      int m=MM[i*width + j];
      // if(m!=0) printf("alpha=%d,beta=%d,m=%d\n",alpha,beta,m);
      for(int u=0; u<m; u++) {
        for(int v=0; v<m-u; v++) {
          // printf("\tu=%d,v=%d\n",u,v);

          //calcul des divergences en alpha,beta
          int k_et = -1;
          for(int k=0; k<=L; k++) {
            discrepancy[k] = D(G, l_polys[k], c, L, u, v, alpha, beta);
            // printf("\t\tk=%d,discrepancy=%d,",k,discrepancy[k]);
            // print_poly(l_polys[k],c,L);
            if(discrepancy[k]!=0) {
              if(k_et==-1) {
                k_et = k;
              } else if(greater_than(l_polys[k_et], l_polys[k], c, L, K)) {
                k_et=k;
              }
            }
          }
          // printf("\tk_et = %d\n",k_et);

          //applications des contraintes en (alpha, beta)
          if(k_et!=-1) { //si la contrainte est vérifiée pour tous les polynomes
            for(int k=0; k<=L; k++) {
              // printf("\t\tk=%d discrepancy=%d\n",k,discrepancy[k]);
              if(discrepancy[k]!=0 && k!=k_et) {
                //Qk=(lambda_ket*Q_k - lambda_k*Q_ket)
                for(int ix=c; ix>=0; ix--) {
                  for(int iy=L; iy>=0; iy--) {
                    l_polys[k][ix][iy] = G->add_table[G->mult_table[discrepancy[k_et]*G->n + l_polys[k][ix][iy]]*G->n + G->mult_table[discrepancy[k]*G->n + l_polys[k_et][ix][iy]]];
                  }
                }
              }
            }
            //Qk_et=(x-alpha)Q_k_et
            for(int ix=c; ix>0; ix--) {
              for(int iy=L; iy>=0; iy--) {
                // l_polys[k_et][ix][iy] = l_polys[k_et][ix-1][iy]; //Q_k*x
                // l_polys[k_et][ix][iy] = G->mult_table[alpha*G->n + l_polys[k_et][ix][iy]];
                l_polys[k_et][ix][iy] = G->add_table[l_polys[k_et][ix-1][iy]*G->n + G->mult_table[alpha*G->n + l_polys[k_et][ix][iy]]];
              }
            }
            for(int iy=0; iy<=L; iy++) l_polys[k_et][0][iy]=G->mult_table[alpha*G->n + l_polys[k_et][0][iy]];
          }

        }
      }
      // for(int cpt=0; cpt<=L && m!=0; cpt++) {
      //   print_poly(l_polys[cpt],c,L);
      // }
    }
  }

  //selection Q
  int k_min=0;
  for(int k=0; k<=L; k++) {
    // print_poly_loc(l_polys[k],c,L);
    if(greater_than(l_polys[k_min], l_polys[k], c, L, K)) k_min=k;
  }

  // print_poly_loc(l_polys[k_min],c,L);
  for(int k=0; k<=L; k++) {
    if(k!=k_min) free_poly(l_polys[k], c);
  }
  u_int8_t** res = l_polys[k_min];
  printf("k_min=%d\n",k_min);
  free(l_polys);

  return res;
}
#include "decoder.h"
#include "string.h"
#include "../generalized_rs/grs.h"

void print_reliability_matrix(double* M, int height, int width) {
  printf("PI =\n");
  for(int i=0; i<height; i++) {
    double s=0;
    for(int j=0; j<width; j++) {
      printf("%.2lf ", M[i*width + j]);
      s+=M[i*width + j];
    }
    printf("| s=%.2lf",s);
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

void print_fword(double* w, int l) {
  printf("w = [%lf",w[0]);
  for(int i=1; i<l; i++) {
    printf(",%lf",w[i]);
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

void print_progressbar(double progress, int total) {
  printf("\rProgress : ");
  int i;
  if(progress>total) progress=total;
  for(i=0; i<progress*20./total; i++) {
    printf("#");
  }
  for(; i<20; i++) {
    printf("-");
  }
  double percent = progress*100.0/(double)total;
  if (percent<0) percent=0.;
  printf(" %.2lf%%",percent);
  fflush(stdout);
}



//msg_len>=2
u_int8_t* decode_soft(double* received, galois* G, int msg_len, int max_mult, double sigma) {
  // print_fword(received, G->n-1);
  u_int8_t* V = get_lambda(G,1);
  u_int8_t* lambda = get_lambda(G,0);
  double* PI = generate_reliabitity_matrix_soft(G, received, sigma);
  // printf("P"); fflush(stdout);
  u_int8_t* M = generate_multiplicity_proportional(G, PI, max_mult, msg_len);
  // printf("M"); fflush(stdout);
  // print_multiplicity_matrix(M, G->n-1, G->n);
  u_int8_t** Q = interpolate(G, M, G->n-1, G->n, msg_len);
  // printf("I"); fflush(stdout);
  //factorization
  int omega = compute_omega(M, G->n-1, G->n, msg_len);
  int L = omega/(msg_len-1);
  int c = cost(M, G->n-1, G->n);
  int n_w;
  int degx_norm, degy_norm;
  u_int8_t** Q_norm = normalize(G, Q, c, L, &degx_norm, &degy_norm);
  u_int8_t** l_w = Factorize(G, Q_norm, degx_norm, degy_norm, 0, msg_len, &n_w);
  // printf("F(%d)",n_w); fflush(stdout);
  if(n_w==0) {
    //décodage non réussi
    free_poly(Q_norm, degx_norm);
    free(l_w);
    u_int8_t* decoded = (u_int8_t*) calloc(msg_len, sizeof(u_int8_t));
    for(int i=0; i<msg_len; i++) {
      decoded[i] = G->n;
    }
    return decoded; //on retourne une valeur impossible
  }
  //selection
  u_int8_t** l_sent = (u_int8_t**) calloc(n_w,sizeof(u_int8_t*));
  for(int i=0;i<n_w;i++) {
    l_sent[i] = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
    evaluate_grs_pol(l_w[i],msg_len,V,lambda,G,l_sent[i]);
  }
  int i_decoded = select_word(l_sent, n_w, G, M);
  u_int8_t* decoded = l_w[i_decoded];
  // printf("S"); fflush(stdout);
  // printf("decoded : ");print_word(l_sent[i_decoded], G->n-1);
  // printf("decoded : ");print_word(decoded, msg_len);
  
  //free
  free_poly(Q_norm, degx_norm);
  for(int i=0; i<n_w; i++) {
    free(l_sent[i]);
    if(i!=i_decoded) free(l_w[i]);
  }
  free(l_w);
  free(l_sent);
  free(M);
  free(PI);
  free(lambda);
  free(V);
  // printf("done\n");

  return decoded;
}
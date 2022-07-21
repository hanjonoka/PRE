#include "reed_solomon/rs.h"
#include "generalized_rs/grs.h"
#include "corps_gallois/gallois.h"
#include "math.h"
#include "stdio.h"

void fprint_poly(FILE* f, u_int8_t* g, int deg){
  for(int i=0;i<deg;i++){
    g[i] ? fprintf(f,"a^%d",g[i]-1) : fprintf(f,"0");
    fprintf(f,"*X^%d + ",i);
  }
  g[deg] ? fprintf(f,"a^%d",g[deg]-1) : fprintf(f,"0");
  fprintf(f,"*X^%d",deg);
  fprintf(f,"\n");
}

void fprint_word(FILE* f, u_int8_t* w, int l){
  for(int i=0; i<l; i++) {
    fprintf(f, "%d", w[i]);
  }
}

void print_word(u_int8_t* w, int l){
  for(int i=0; i<l; i++) {
    printf("%d", w[i]);
  }
  printf("\n");
}

int main() {
  u_int8_t P[3] =  {1,1,1};
  galois* G = generate_galois(P, 2);
  int msg_len = 2;

  u_int8_t** msg_tab = (u_int8_t**) calloc(pow(G->n,msg_len),sizeof(u_int8_t*));

  //init RS
  u_int8_t** g_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  int deg_g = G->n-1-msg_len;
  u_int8_t* racines = calloc(deg_g,sizeof(u_int8_t));
  for(int i=0;i<G->n-1;i++) {
    for(int j=0;j<deg_g;j++) racines[j] = ((i+j)%(G->n-1)+1);
    g_tab[i] = get_g(G,racines,deg_g);
  }

  //init GRS
  u_int8_t* V = get_lambda(G,1);
  u_int8_t** lambda_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int i=0; i<G->n-1; i++) {
    lambda_tab[i] = get_lambda(G, i);
  }

  FILE* rs4_f = fopen("rs4.csv","w");
  FILE* grs4_f = fopen("grs4.csv","w");
  fprintf(rs4_f,"sep=,\n");
  fprintf(grs4_f,"sep=,\n");

  u_int8_t** rs_enc_tab = (u_int8_t**)calloc(pow(G->n,msg_len),sizeof(u_int8_t*));
  u_int8_t** grs_enc_tab = (u_int8_t**)calloc(pow(G->n,msg_len),sizeof(u_int8_t*));

  // printf("init done\n");
  for(u_int16_t i=0;i<G->n-1;i++) {
    fprint_poly(rs4_f,g_tab[i],deg_g);
    fprintf(grs4_f,"lambda=");fprint_word(grs4_f,lambda_tab[i], G->n-1);fprintf(grs4_f,"\n");
    // printf("i=%d\n",i);
    for(u_int8_t j=0; j<pow(G->n,msg_len); j++) {
      int n=j;
      msg_tab[j] = (u_int8_t*) calloc(msg_len, sizeof(u_int8_t));
      rs_enc_tab[j] = realloc(rs_enc_tab[j], (G->n-1)*sizeof(u_int8_t));
      grs_enc_tab[j] = realloc(grs_enc_tab[j], (G->n-1)*sizeof(u_int8_t));
      // printf("init msg\n");
      for(int k = 0; k<msg_len; k++) {
        msg_tab[j][k] = n%G->n;
        n=n/G->n;
      }
      // print_word(msg_tab[j], msg_len);
      evaluate_rs_pol(msg_tab[j],msg_len,g_tab[i],G,rs_enc_tab[j]);
      evaluate_grs_pol(msg_tab[j], msg_len, V, lambda_tab[i], G, grs_enc_tab[j]);
      // print_word(rs_enc_tab[j],G->n-1);
    }

    for(int j=0;j<pow(G->n,msg_len);j++){
      fprint_word(rs4_f, msg_tab[j], msg_len);
      fprintf(rs4_f,",");
      fprint_word(rs4_f, rs_enc_tab[j], G->n-1);
      fprintf(rs4_f,"\n");

      fprint_word(grs4_f, msg_tab[j], msg_len);
      fprintf(grs4_f,",");
      fprint_word(grs4_f, grs_enc_tab[j], G->n-1);
      fprintf(grs4_f,"\n");
    }
    // printf("--------\n");

    fprintf(rs4_f,"\n");
    fprintf(grs4_f,"\n");
    fflush(rs4_f);
    fflush(grs4_f);
  }

  fclose(rs4_f);
  fclose(grs4_f);

  return 0;
}
#include "reed_solomon/rs.h"
#include "generalized_rs/grs.h"
#include "corps_gallois/gallois.h"
#include "math.h"
#include "stdio.h"

int poids_mot(u_int8_t* m, int len) {
  int p = 0;
  for(int i=0;i<len;i++) {
    if(m[i]) p++;
  }
  return p;
}

int poids_bits_mot(u_int8_t* m, int len, galois* G) {
  int p = 0;
  for(int i=0;i<len;i++) p+=G->poids[m[i]];
  return p;
}

int main() {
  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);

  //init RS
  u_int8_t** g_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int i=0;i<G->n-1;i++) {
    u_int8_t racines[2] = {i+1,(i+1)%(G->n-1)+1};
    g_tab[i] = get_g(G,racines,2);
  }

  //init GRS
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t** lambda_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int i=0;i<G->n-1;i++) {
    lambda_tab[i] = get_lambda(G,i);
  }

  FILE* wh_sym = fopen("wh_sym.csv","w");
  fprintf(wh_sym,"sep=,\n");
  fprintf(wh_sym,"poids");
  for(int i=0;i<G->n-1;i++)  fprintf(wh_sym,",(x-a^%d)(x-a^%d)",i,(i)%(G->n-1)+1);
  for(int i=0;i<G->n-1;i++)  fprintf(wh_sym,",a^(%d)",i);
  fprintf(wh_sym,"\n");

  FILE* wh_bit = fopen("wh_bit.csv","w");
  fprintf(wh_bit,"sep=,\n");
  fprintf(wh_bit,"poids");
  for(int i=0;i<G->n-1;i++)  fprintf(wh_bit,",(x-a^%d)(x-a^%d)",i,(i+1)%(G->n-1));
  for(int i=0;i<G->n-1;i++)  fprintf(wh_bit,",a^(%d)",i);
  fprintf(wh_bit,"\n");

  // u_int16_t* poids_grs1 = (u_int16_t*)calloc(8,sizeof(u_int16_t));
  // u_int16_t* poids_grs2 = (u_int16_t*)calloc(8,sizeof(u_int16_t));
  // u_int16_t* poids_grs3 = (u_int16_t*)calloc(8,sizeof(u_int16_t));

  // u_int16_t* poids_bits_grs1 = (u_int16_t*)calloc(7*3+1,sizeof(u_int16_t));
  // u_int16_t* poids_bits_grs2 = (u_int16_t*)calloc(7*3+1,sizeof(u_int16_t));
  // u_int16_t* poids_bits_grs3 = (u_int16_t*)calloc(7*3+1,sizeof(u_int16_t));

  u_int16_t** poids_rs = (u_int16_t**)calloc(G->n-1,sizeof(u_int16_t*));
  u_int16_t** poids_bits_rs = (u_int16_t**)calloc(G->n-1,sizeof(u_int16_t*));
  for(int i=0;i<G->n-1;i++) {
    poids_rs[i] = (u_int16_t*)calloc(G->n,sizeof(u_int16_t));
    poids_bits_rs[i] = (u_int16_t*)calloc((G->n-1)*(G->deg_P)+1,sizeof(u_int16_t));
  }

  u_int16_t** poids_grs = (u_int16_t**)calloc(G->n-1,sizeof(u_int16_t*));
  u_int16_t** poids_bits_grs = (u_int16_t**)calloc(G->n-1,sizeof(u_int16_t*));
  for(int i=0;i<G->n-1;i++) {
    poids_grs[i] = (u_int16_t*)calloc(G->n,sizeof(u_int16_t));
    poids_bits_grs[i] = (u_int16_t*)calloc((G->n-1)*(G->deg_P)+1,sizeof(u_int16_t));
  }

  int msg_len = 5;
  u_int8_t* msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t));
  u_int8_t* coded = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
  for(u_int16_t i=0;i<pow(G->n,msg_len);i++) {
    int n=i;
    for(int j=0;j<G->n-G->deg_P;j++) {
      int x = (n % ((int) pow(G->n,(j+1))));
      n=n-x;
      msg[j]=x/pow(G->n,j);
    }

    for(int i=0;i<G->n-1;i++) {
      evaluate_rs_pol(msg,msg_len,g_tab[i],G,coded);
      poids_rs[i][poids_mot(coded,G->n-1)]++;
      poids_bits_rs[i][poids_bits_mot(coded,G->n-1,G)]++;

      evaluate_grs_pol(msg,msg_len,V,lambda_tab[i],G,coded);
      poids_grs[i][poids_mot(coded,G->n-1)]++;
      poids_bits_grs[i][poids_bits_mot(coded,G->n-1,G)]++;
    }

    // evaluate_grs_pol(msg,5,V,lambda1,G,coded);
    // poids_grs1[poids_mot(coded,7)]++;
    // poids_bits_grs1[poids_bits_mot(coded,7,G)]++;
    // evaluate_grs_pol(msg,5,V,lambda2,G,coded);
    // poids_grs2[poids_mot(coded,7)]++;
    // poids_bits_grs2[poids_bits_mot(coded,7,G)]++;
    // evaluate_grs_pol(msg,5,V,lambda3,G,coded);
    // poids_grs3[poids_mot(coded,7)]++;
    // poids_bits_grs3[poids_bits_mot(coded,7,G)]++;

  }

  for(int i=0;i<G->n;i++) {
    fprintf(wh_sym,"%d",i);
    for(int j=0;j<G->n-1;j++) {
      fprintf(wh_sym,",%d",poids_rs[j][i]);
    }
    for(int j=0;j<G->n-1;j++) {
      fprintf(wh_sym,",%d",poids_grs[j][i]);
    }
    fprintf(wh_sym,"\n");
    // fprintf(wh_sym,"%d,%d,%d,%d,%d,%d,%d\n",i,poids_rs1[i],poids_rs2[i],poids_rs3[i],poids_grs1[i],poids_grs2[i],poids_grs3[i]);
  }

  for(int i=0;i<(G->n-1)*(G->deg_P)+1;i++) {
    fprintf(wh_bit,"%d",i);
    for(int j=0;j<G->n-1;j++) {
      fprintf(wh_bit,",%d",poids_bits_rs[j][i]);
    }
    for(int j=0;j<G->n-1;j++) {
      fprintf(wh_bit,",%d",poids_bits_grs[j][i]);
    }
    fprintf(wh_bit,"\n");
    // fprintf(wh_sym,"%d,%d,%d,%d,%d,%d,%d\n",i,poids_rs1[i],poids_rs2[i],poids_rs3[i],poids_grs1[i],poids_grs2[i],poids_grs3[i]);
  }


  fclose(wh_sym);
  fclose(wh_bit);

  return 0;
}
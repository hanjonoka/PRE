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
  u_int8_t P[5] =  {1,1,0,0,1};
  galois* G = generate_galois(P, 4);

  //init RS
  u_int8_t racines1[6] = {2,3,4,5,6,7};
  u_int8_t* g1 = get_g(G, racines1, 6);
  u_int8_t racines2[6] = {3,4,5,6,7,8};
  u_int8_t* g2 = get_g(G, racines2, 6);
  u_int8_t racines3[6] = {4,5,6,7,8,9};
  u_int8_t* g3 = get_g(G, racines3, 6);

  //init GRS
  u_int8_t V[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t* lambda1 = get_lambda(G,15,0); //coeficients de normalisation.
  u_int8_t* lambda2 = get_lambda(G,7,1); //coeficients de normalisation.
  u_int8_t* lambda3 = get_lambda(G,7,2); //coeficients de normalisation.

  FILE* wh_sym = fopen("wh_sym16.csv","w");
  fprintf(wh_sym,"sep=,\n");
  fprintf(wh_sym,"poids,{2,3,4,5,6,7},{3,4,5,6,7,8},{4,5,6,7,8,9},a1,a2,a3\n");

  FILE* wh_bit = fopen("wh_bit16.csv","w");
  fprintf(wh_bit,"sep=,\n");
  fprintf(wh_bit,"poids,{2,3,4,5,6,7},{3,4,5,6,7,8},{4,5,6,7,8,9},a1,a2,a3\n");

  int code_len = 15;
  int msg_len = 9;

  u_int16_t* poids_rs1 = (u_int16_t*)calloc(code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_rs2 = (u_int16_t*)calloc(code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_rs3 = (u_int16_t*)calloc(code_len+1,sizeof(u_int16_t));

  u_int16_t* poids_grs1 = (u_int16_t*)calloc(code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_grs2 = (u_int16_t*)calloc(code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_grs3 = (u_int16_t*)calloc(code_len+1,sizeof(u_int16_t));

  u_int16_t* poids_bits_rs1 = (u_int16_t*)calloc(G->deg_P*code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_bits_rs2 = (u_int16_t*)calloc(G->deg_P*code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_bits_rs3 = (u_int16_t*)calloc(G->deg_P*code_len+1,sizeof(u_int16_t));

  u_int16_t* poids_bits_grs1 = (u_int16_t*)calloc(G->deg_P*code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_bits_grs2 = (u_int16_t*)calloc(G->deg_P*code_len+1,sizeof(u_int16_t));
  u_int16_t* poids_bits_grs3 = (u_int16_t*)calloc(G->deg_P*code_len+1,sizeof(u_int16_t));

  u_int8_t* msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t));
  u_int8_t* coded = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t));

  for(u_int64_t i=0;i<pow(G->n,msg_len);i++) {
    if(i%1000000==0) printf("%ld\n",i);
    int n=i;
    for(int j=0;j<msg_len;j++) {
      int x = (n % ((int) pow(G->n,(j+1))));
      n=n-x;
      msg[j]=x/pow(G->n,j);
    }

    evaluate_rs_pol(msg,msg_len,g1,G,coded);
    poids_rs1[poids_mot(coded,7)]++;
    poids_bits_rs1[poids_bits_mot(coded,7,G)]++;
    // evaluate_rs_pol(msg,msg_len,g2,G,coded);
    // poids_rs2[poids_mot(coded,7)]++;
    // poids_bits_rs2[poids_bits_mot(coded,7,G)]++;
    // evaluate_rs_pol(msg,msg_len,g3,G,coded);
    // poids_rs3[poids_mot(coded,7)]++;
    // poids_bits_rs3[poids_bits_mot(coded,7,G)]++;

    // evaluate_grs_pol(msg,msg_len,V,lambda1,G,coded);
    // poids_grs1[poids_mot(coded,7)]++;
    // poids_bits_grs1[poids_bits_mot(coded,7,G)]++;
    // evaluate_grs_pol(msg,msg_len,V,lambda2,G,coded);
    // poids_grs2[poids_mot(coded,7)]++;
    // poids_bits_grs2[poids_bits_mot(coded,7,G)]++;
    // evaluate_grs_pol(msg,msg_len,V,lambda3,G,coded);
    // poids_grs3[poids_mot(coded,7)]++;
    poids_bits_grs3[poids_bits_mot(coded,7,G)]++;
  }

  for(int i=0;i<G->n;i++) {
    fprintf(wh_sym,"%d,%d,%d,%d,%d,%d,%d\n",i,poids_rs1[i],poids_rs2[i],poids_rs3[i],poids_grs1[i],poids_grs2[i],poids_grs3[i]);
  }

  for(int i=0;i<G->deg_P*code_len+1;i++) {
    fprintf(wh_bit,"%d,%d,%d,%d,%d,%d,%d\n",i,poids_bits_rs1[i],poids_bits_rs2[i],poids_bits_rs3[i],poids_bits_grs1[i],poids_bits_grs2[i],poids_bits_grs3[i]);
  }

  fclose(wh_sym);
  fclose(wh_bit);

  return 0;
}
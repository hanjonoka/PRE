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
  u_int8_t racines1[2] = {1,2};
  u_int8_t* g1 = get_g(G, racines1, 2);
  u_int8_t racines2[2] = {2,3};
  u_int8_t* g2 = get_g(G, racines2, 2);
  u_int8_t racines3[2] = {3,4};
  u_int8_t* g3 = get_g(G, racines3, 2);

  //init GRS
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t* lambda1 = get_lambda(G,7,1); //coeficients de normalisation.
  u_int8_t* lambda2 = get_lambda(G,7,0); //coeficients de normalisation.
  u_int8_t* lambda3 = get_lambda(G,7,5); //coeficients de normalisation.

  FILE* wh_sym = fopen("wh_sym.csv","w");
  fprintf(wh_sym,"sep=,\n");
  fprintf(wh_sym,"poids,(x-1)(x-a),(x-a)(x-a2),(x-a2)(x-a4),a1,a0,a5\n");

  FILE* wh_bit = fopen("wh_bit.csv","w");
  fprintf(wh_bit,"sep=,\n");
  fprintf(wh_bit,"poids,(x-1)(x-a),(x-a)(x-a2),(x-a2)(x-a4),a1,a0,a5\n");

  u_int16_t* poids_rs1 = (u_int16_t*)calloc(8,sizeof(u_int16_t));
  u_int16_t* poids_rs2 = (u_int16_t*)calloc(8,sizeof(u_int16_t));
  u_int16_t* poids_rs3 = (u_int16_t*)calloc(8,sizeof(u_int16_t));

  u_int16_t* poids_grs1 = (u_int16_t*)calloc(8,sizeof(u_int16_t));
  u_int16_t* poids_grs2 = (u_int16_t*)calloc(8,sizeof(u_int16_t));
  u_int16_t* poids_grs3 = (u_int16_t*)calloc(8,sizeof(u_int16_t));

  u_int16_t* poids_bits_rs1 = (u_int16_t*)calloc(22,sizeof(u_int16_t));
  u_int16_t* poids_bits_rs2 = (u_int16_t*)calloc(22,sizeof(u_int16_t));
  u_int16_t* poids_bits_rs3 = (u_int16_t*)calloc(22,sizeof(u_int16_t));

  u_int16_t* poids_bits_grs1 = (u_int16_t*)calloc(7*3+1,sizeof(u_int16_t));
  u_int16_t* poids_bits_grs2 = (u_int16_t*)calloc(7*3+1,sizeof(u_int16_t));
  u_int16_t* poids_bits_grs3 = (u_int16_t*)calloc(7*3+1,sizeof(u_int16_t));

  u_int8_t* msg = (u_int8_t*) calloc(5,sizeof(u_int8_t));
  u_int8_t* coded = (u_int8_t*) calloc(7,sizeof(u_int8_t));
  for(u_int16_t i=0;i<pow(8,5);i++) {
    int n=i;
    for(int j=0;j<5;j++) {
      int x = (n % ((int) pow(8,(j+1))));
      n=n-x;
      msg[j]=x/pow(8,j);
    }

    evaluate_rs_pol(msg,5,g1,G,coded);
    poids_rs1[poids_mot(coded,7)]++;
    poids_bits_rs1[poids_bits_mot(coded,7,G)]++;
    evaluate_rs_pol(msg,5,g2,G,coded);
    poids_rs2[poids_mot(coded,7)]++;
    poids_bits_rs2[poids_bits_mot(coded,7,G)]++;
    evaluate_rs_pol(msg,5,g3,G,coded);
    poids_rs3[poids_mot(coded,7)]++;
    poids_bits_rs3[poids_bits_mot(coded,7,G)]++;

    evaluate_grs_pol(msg,5,V,lambda1,G,coded);
    poids_grs1[poids_mot(coded,7)]++;
    poids_bits_grs1[poids_bits_mot(coded,7,G)]++;
    evaluate_grs_pol(msg,5,V,lambda2,G,coded);
    poids_grs2[poids_mot(coded,7)]++;
    poids_bits_grs2[poids_bits_mot(coded,7,G)]++;
    evaluate_grs_pol(msg,5,V,lambda3,G,coded);
    poids_grs3[poids_mot(coded,7)]++;
    poids_bits_grs3[poids_bits_mot(coded,7,G)]++;

  }

  for(int i=0;i<8;i++) {
    fprintf(wh_sym,"%d,%d,%d,%d,%d,%d,%d\n",i,poids_rs1[i],poids_rs2[i],poids_rs3[i],poids_grs1[i],poids_grs2[i],poids_grs3[i]);
  }

  for(int i=0;i<7*3+1;i++) {
    fprintf(wh_bit,"%d,%d,%d,%d,%d,%d,%d\n",i,poids_bits_rs1[i],poids_bits_rs2[i],poids_bits_rs3[i],poids_bits_grs1[i],poids_bits_grs2[i],poids_bits_grs3[i]);
  }

  fclose(wh_sym);
  fclose(wh_bit);

  return 0;
}
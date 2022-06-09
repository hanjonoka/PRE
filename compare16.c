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
  int msg_len = 7; // CHANGER NOM FICHIER AUSSI

  //init RS
  u_int8_t** g_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  int deg_g = G->n-1-msg_len;
  u_int8_t* racines = calloc(deg_g,sizeof(u_int8_t));
  for(int i=0;i<G->n-1;i++) {
    // u_int8_t racines_old[10] = {i+1,(i+1)%(G->n-1)+1,(i+2)%(G->n-1)+1,(i+3)%(G->n-1)+1,(i+4)%(G->n-1)+1,(i+5)%(G->n-1)+1,(i+6)%(G->n-1)+1,(i+7)%(G->n-1)+1,(i+8)%(G->n-1)+1,(i+9)%(G->n-1)+1};
    for(int j=0;j<deg_g;j++) racines[j] = ((i+j)%(G->n-1)+1);
    g_tab[i] = get_g(G,racines,deg_g);
    printf("%d : ",i);print_poly(g_tab[i],deg_g);
  }
  free(racines);

  //init GRS
  // u_int8_t V[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t* V = get_lambda(G,1);
  u_int8_t** lambda_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int i=0;i<G->n-1;i++) {
    lambda_tab[i] = get_lambda(G,i);
    printf("%d : ",i);print_lambda(lambda_tab[i],G->n-1);
  }

  FILE* wh_sym = fopen("wh_sym16_7-15.csv","w");
  fprintf(wh_sym,"sep=,\n");
  fprintf(wh_sym,"poids");
  for(int i=0;i<G->n-1;i++)  fprintf(wh_sym,",(x-a^%d):(x-a^%d)",i,(i+9)%(G->n-1));
  for(int i=0;i<G->n-1;i++)  fprintf(wh_sym,",a^(%d)",i);
  fprintf(wh_sym,"\n");

  FILE* wh_bit = fopen("wh_bit16_7-15.csv","w");
  fprintf(wh_bit,"sep=,\n");
  fprintf(wh_bit,"poids");
  for(int i=0;i<G->n-1;i++)  fprintf(wh_bit,",(x-a^%d):(x-a^%d)",i,(i+9)%(G->n-1));
  for(int i=0;i<G->n-1;i++)  fprintf(wh_bit,",a^(%d)",i);
  fprintf(wh_bit,"\n");


  u_int32_t** poids_rs = (u_int32_t**)calloc(G->n-1,sizeof(u_int32_t*));
  u_int32_t** poids_bits_rs = (u_int32_t**)calloc(G->n-1,sizeof(u_int32_t*));
  for(int i=0;i<G->n-1;i++) {
    poids_rs[i] = (u_int32_t*)calloc(G->n,sizeof(u_int32_t));
    poids_bits_rs[i] = (u_int32_t*)calloc((G->n-1)*(G->deg_P)+1,sizeof(u_int32_t));
  }

  u_int32_t** poids_grs = (u_int32_t**)calloc(G->n-1,sizeof(u_int32_t*));
  u_int32_t** poids_bits_grs = (u_int32_t**)calloc(G->n-1,sizeof(u_int32_t*));
  for(int i=0;i<G->n-1;i++) {
    poids_grs[i] = (u_int32_t*)calloc(G->n,sizeof(u_int32_t));
    poids_bits_grs[i] = (u_int32_t*)calloc((G->n-1)*(G->deg_P)+1,sizeof(u_int32_t));
  }

  u_int8_t* msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t));
  u_int8_t* coded = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
  printf("imax = %d\n",(int)pow(G->n,msg_len));
  for(u_int32_t i=0;i<pow(G->n,msg_len);i++) {
    if(i%10000==0) printf("%d\n",i);
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
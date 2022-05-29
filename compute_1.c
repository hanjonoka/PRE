#include "reed_solomon/rs.h"
#include "generalized_rs/grs.h"
#include "corps_gallois/gallois.h"
#include "math.h"
#include "stdio.h"

int main() {
  gallois* G = generate_gallois_8();

  //init RS
  u_int8_t racines[2] = {2,3};
  u_int8_t* g = get_g(G, racines, 2);

  //init GRS
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda[7] = {1,1,1,1,1,1,1}; //coeficients de normalisation.

  FILE* rs_f = fopen("rs.csv","w");
  FILE* grs_f = fopen("grs.csv","w");
  fprintf(rs_f,"sep=,\n");
  fprintf(grs_f,"sep=,\n");

  u_int8_t** msg = (u_int8_t**)calloc(5,sizeof(u_int8_t*));
  for (int i=0;i<5;i++) {
    msg[i] = (u_int8_t*)calloc(3,sizeof(u_int8_t));
  }

  u_int8_t** rs_enc = (u_int8_t**)calloc(7,sizeof(u_int8_t*));
  u_int8_t** grs_enc = (u_int8_t**)calloc(7,sizeof(u_int8_t*));

  for(u_int16_t i=0;i<pow(8,5);i++) {
    if(i%1000==0) printf("i = %d\n",i);
    for(int j=0;j<5;j++){
      for(int k=0;k<3;k++){
        int n = ((4-j)*3 + (2-k));
        msg[j][k] = (i>>n) % 2;
        fprintf(rs_f,"%d", msg[j][k]);
        fprintf(grs_f,"%d", msg[j][k]);
      }
    }
    fprintf(rs_f,",");
    fprintf(grs_f,",");

    evaluate_grs(msg,5,V,lambda,G,grs_enc);
    evaluate_rs(msg,5,g,G,rs_enc);

    for(int j=0;j<7;j++){
      for(int k=0;k<3;k++){
        fprintf(rs_f,"%d", rs_enc[j][k]);
        fprintf(grs_f,"%d", grs_enc[j][k]);
      }
    }

    fprintf(rs_f,"\n");
    fprintf(grs_f,"\n");
  }

  fclose(rs_f);
  fclose(grs_f);

  return 0;
}
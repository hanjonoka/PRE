#include "reed_solomon/rs.h"
#include "generalized_rs/grs.h"
#include "corps_gallois/gallois.h"
#include "math.h"
#include "stdio.h"

int main() {
  //----8 Elements---------------
  u_int8_t P8[4] =  {1,0,1,1};
  galois* G8 = generate_galois(P8, 3);

  //init RS
  u_int8_t racines8[2] = {2,3};
  u_int8_t* g8 = get_g(G8, racines8, 2);

  //init GRS
  u_int8_t V8[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda8[7] = {1,1,1,1,1,1,1}; //coeficients de normalisation.

  FILE* rs8_f = fopen("rs8.csv","w");
  FILE* grs8_f = fopen("grs8.csv","w");
  fprintf(rs8_f,"sep=,\n");
  fprintf(grs8_f,"sep=,\n");

  u_int8_t** msg8 = (u_int8_t**)calloc(5,sizeof(u_int8_t*));
  for (int i=0;i<5;i++) {
    msg8[i] = (u_int8_t*)calloc(3,sizeof(u_int8_t));
  }

  u_int8_t** rs_enc = (u_int8_t**)calloc(7,sizeof(u_int8_t*));
  u_int8_t** grs_enc = (u_int8_t**)calloc(7,sizeof(u_int8_t*));

  for(u_int16_t i=0;i<pow(8,5);i++) {
    if(i%1000==0) printf("i = %d\n",i);
    for(int j=0;j<5;j++){
      for(int k=0;k<3;k++){
        int n = ((4-j)*3 + (2-k));
        msg8[j][k] = (i>>n) % 2;
        fprintf(rs8_f,"%d", msg8[j][k]);
        fprintf(grs8_f,"%d", msg8[j][k]);
      }
    }
    fprintf(rs8_f,",");
    fprintf(grs8_f,",");

    evaluate_grs(msg8,5,V8,lambda8,G8,grs_enc);
    evaluate_rs(msg8,5,g8,G8,rs_enc);

    for(int j=0;j<7;j++){
      for(int k=0;k<3;k++){
        fprintf(rs8_f,"%d", rs_enc[j][k]);
        fprintf(grs8_f,"%d", grs_enc[j][k]);
      }
    }

    fprintf(rs8_f,"\n");
    fprintf(grs8_f,"\n");
  }

  fclose(rs8_f);
  fclose(grs8_f);

  //--16 elements-------------------------------------------
  u_int8_t P16[5] =  {1,1,0,0,1};
  galois* G16 = generate_galois(P16, 4);

  //init RS
  u_int8_t racines16[6] = {2,3,4,5,6,7};
  u_int8_t* g16 = get_g(G16, racines16, 6);

  //init GRS
  u_int8_t V16[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda16[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; //coeficients de normalisation.

  FILE* rs16_f = fopen("rs16.csv","w");
  FILE* grs16_f = fopen("grs16.csv","w");
  fprintf(rs8_f,"sep=,\n");
  fprintf(grs8_f,"sep=,\n");

  u_int8_t** msg16 = (u_int8_t**)calloc(9,sizeof(u_int8_t*));
  for (int i=0;i<9;i++) {
    msg16[i] = (u_int8_t*)calloc(4,sizeof(u_int8_t));
  }

  u_int8_t** rs16_enc = (u_int8_t**)calloc(15,sizeof(u_int8_t*));
  u_int8_t** grs16_enc = (u_int8_t**)calloc(15,sizeof(u_int8_t*));

  for(u_int64_t i=0;i<100;i++) {
    printf("i = %ld\n",i);
    for(int j=0;j<9;j++){
      for(int k=0;k<4;k++){
        int n = ((8-j)*4 + (3-k));
        msg16[j][k] = (i>>n) % 2;
        fprintf(rs16_f,"%d", msg16[j][k]);
        fprintf(grs16_f,"%d", msg16[j][k]);
      }
    }
    fprintf(rs16_f,",");
    fprintf(grs16_f,",");

    evaluate_grs(msg16,9,V16,lambda16,G16,grs16_enc);
    evaluate_rs(msg16,9,g16,G16,rs16_enc);

    for(int j=0;j<15;j++){
      for(int k=0;k<4;k++){
        fprintf(rs16_f,"%d", rs16_enc[j][k]);
        fprintf(grs16_f,"%d", grs16_enc[j][k]);
      }
    }

    fprintf(rs16_f,"\n");
    fprintf(grs16_f,"\n");
  }
  fclose(rs16_f);
  fclose(grs16_f);


  return 0;
}
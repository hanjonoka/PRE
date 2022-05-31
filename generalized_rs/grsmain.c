#include "grs.h"


int main(){
  printf("----> Generalized Reed Solomon 8\n");
  u_int8_t msg_h[5][3] = {{0,0,0},{0,1,0},{0,0,0},{0,0,0},{1,0,1}};
  u_int8_t** msg = (u_int8_t**) calloc(5,sizeof(u_int8_t*));
  for(int i=0;i<5;i++) {
    msg[i] = (u_int8_t*) calloc(3,sizeof(u_int8_t));
    for (int j=0;j<3;j++){
      msg[i][j] =  msg_h[i][j];
    }
  }

  printf("message : [");
  for(int i=0;i<5;i++){
    printf("[%d,%d,%d]",msg[i][0],msg[i][1],msg[i][2]);
  }
  printf("]\n");


  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda[7] = {2,2,2,2,2,2,2}; //coeficients de normalisation.
  u_int8_t** send_msg = (u_int8_t**)malloc(7*sizeof(u_int8_t*));
  evaluate_grs(msg,5,V,lambda,G,send_msg);

  printf("encoded : [");
  for(int i=0;i<7;i++){
    printf("[%d,%d,%d]",send_msg[i][0],send_msg[i][1],send_msg[i][2]);
  }
  printf("\n");

  printf("----> Generalized Reed Solomon 16\n");
  u_int8_t msg2_h[9][4] = {{1,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
  u_int8_t** msg2 = (u_int8_t**)calloc(9,sizeof(u_int8_t*));
  for(int i=0;i<9;i++) {
    msg2[i] = (u_int8_t*)calloc(4,sizeof(u_int8_t));
    for (int j=0;j<4;j++){
      msg2[i][j] =  msg2_h[i][j];
    }
  }
  printf("message : [");
  for(int i=0;i<9;i++){
    printf("[%d,%d,%d,%d]",msg2[i][0],msg2[i][1],msg2[i][2],msg2[i][3]);
  }
  printf("]\n");

  u_int8_t P2[5] = {1,1,0,0,1};
  galois* G2 = generate_galois(P2,4);
  u_int8_t V2[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda2[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; //coeficients de normalisation.
  u_int8_t** send_msg2 = (u_int8_t**)malloc(15*sizeof(u_int8_t*));
  evaluate_grs(msg2,9,V2,lambda2,G2,send_msg2);

  printf("encoded : [");
  for(int i=0;i<G2->n-1;i++){
    printf("[%d,%d,%d,%d]",send_msg2[i][0],send_msg2[i][1],send_msg2[i][2],send_msg2[i][3]);
  }
  printf("]\n");

  return 0;
}
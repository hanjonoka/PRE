#include "rs.h"

int main(){

  u_int8_t msg_h[5][3] = {{0,0,0},{0,1,0},{0,0,0},{0,0,0},{1,0,1}};
  u_int8_t** msg = (u_int8_t**)calloc(5,sizeof(u_int8_t*));
  for(int i=0;i<5;i++) {
    msg[i] = (u_int8_t*)calloc(3,sizeof(u_int8_t));
    for (int j=0;j<3;j++){
      msg[i][j] =  msg_h[i][j];
    }
  }

  printf("message : [");
  for(int i=0;i<5;i++){
    printf("[%d,%d,%d]",msg[i][0],msg[i][1],msg[i][2]);
  }
  printf("]\n");

  printf("----> Reed Solomon\n");
  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  u_int8_t racines[2] = {2,3};
  u_int8_t* g = get_g(G, racines, 2);
  print_poly(g,2);
  u_int8_t** send_msg = (u_int8_t**) malloc(7*sizeof(u_int8_t*));
  evaluate_rs(msg,5,g,G,send_msg);
  printf("encoded : [");
  for(int i=0;i<G->n-1;i++){
    printf("[%d,%d,%d]",send_msg[i][0],send_msg[i][1],send_msg[i][2]);
  }
  printf("]\n");

//--------------------------------------------------------------
printf("----------------------------------------------------\n");

  u_int8_t msg2_h[9][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}};
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
  u_int8_t racines2[6] = {2,3,4,5,6,7};
  u_int8_t* g2 = get_g(G2,racines2,6);
  printf("g : "); print_poly(g2,6);
  u_int8_t** send_msg2 = (u_int8_t**) malloc((G2->n-1)*sizeof(u_int8_t*));

  evaluate_rs(msg2,9,g2,G2,send_msg2);
  printf("encoded : [");
  for(int i=0;i<G2->n-1;i++){
    printf("[%d,%d,%d,%d]",send_msg2[i][0],send_msg2[i][1],send_msg2[i][2],send_msg2[i][3]);
  }
  printf("]\n");



  return 0;

}
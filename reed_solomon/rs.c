#include "rs.h"

//tableau de formes polaires de G8 (u_int8_t)
u_int8_t* get_g(galois* G, u_int8_t* elts, int deg){
  u_int8_t* g = (u_int8_t*) calloc(deg+1,sizeof(u_int8_t));
  g[0]=1;
  u_int8_t* tmp = (u_int8_t*) calloc(deg+1,sizeof(u_int8_t));
  for(int i=0;i<deg;i++){
    for(int j=0;j<=deg;j++){
      tmp[j] = G->add_table[G->mult_table[g[j]*G->n + (elts[i])]*G->n + (j>0 ? g[j-1] : 0)];
    }
    for(int j=0;j<=deg;j++) g[j] = tmp[j];
  }
  free(tmp);
  return g;
}

void print_poly(u_int8_t* g, int deg){
  for(int i=0;i<deg;i++){
    g[i] ? printf("a^%d",g[i]-1) : printf("0");
    printf("*X^%d + ",i);
  }
  g[deg] ? printf("a^%d",g[deg]-1) : printf("0");
  printf("*X^%d",deg);
  printf("\n");
}


void evaluate_rs_old(u_int8_t** msg, u_int8_t* g, galois* G, u_int8_t** encoded) {
  u_int8_t* c_msg = (u_int8_t*)calloc(7,sizeof(u_int8_t));
  u_int8_t* r = (u_int8_t*)calloc(7,sizeof(u_int8_t));
  for(int i=2;i<7;i++) {
    c_msg[i]=pol_of_cart_8(G->elements,msg[i-2]);
    r[i]=c_msg[i];
  }

  //calcul du reste mod g
  print_poly(r,6);
  for(int i=6; i>=2; i--){
    for(int j=i-2; j<=i; j++){
      u_int8_t m = (G->mult_table[r[i]*8 + g[j-i+2]]);
      // printf("i:%d j:%d r[j]:%d j-i+2:%d m:%d\n",i,j,r[j],j-i+2,m);
      r[j] = G->add_table[r[j]*8 + m];
    }
  }

  //adding r to msg
  c_msg[0]=r[0];
  c_msg[1]=r[1];

  for(int i=0;i<7;i++){
    int c = c_msg[i];
    encoded[i]=cart_of_pol(G, c);
  }

  free(r);
  free(c_msg);
}

void evaluate_rs(u_int8_t** msg, int msg_len, u_int8_t* g, galois* G, u_int8_t** encoded) {
  // int encoded_len = G->n-1;
  // int deg_g = encoded_len-msg_len;
  // u_int8_t* c_msg = (u_int8_t*)calloc(encoded_len, sizeof(u_int8_t));
  // u_int8_t* r = (u_int8_t*)calloc(encoded_len, sizeof(u_int8_t));
  // for(int i=deg_g; i<encoded_len; i++) {
  //   c_msg[i]=pol_of_cart(G,msg[i-deg_g]);
  //   r[i] = c_msg[i];
  // }

  // //calcul du reste mod g
  // for(int i=encoded_len-1; i>=deg_g; i--) {
  //   for(int j=i-deg_g; j<=i; j++) {
  //     u_int8_t m = (G->mult_table[r[i]*G->n + g[j-i+deg_g]]);
  //     // printf("i:%d j:%d r[j]:%d j-i+2:%d m:%d\n",i,j,r[j],j-i+2,m);
  //     r[j] = G->add_table[r[j]*G->n + m];
  //   }
  // }

  // //adding r to msg
  // for(int i=0; i<deg_g; i++) c_msg[i] = r[i];

  // for(int i=0; i<encoded_len;i++) {
  //   encoded[i] = cart_of_pol(G, c_msg[i]);
  // }
  // // print_poly(c_msg,encoded_len-1);

  // free(r);
  // free(c_msg);
  u_int8_t* pol_msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t));
  for(int i=0;i<msg_len;i++) pol_msg[i]=pol_of_cart(G,msg[i]);
  u_int8_t* pol_encoded = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
  evaluate_rs_pol(pol_msg,msg_len,g,G,pol_encoded);
  for(int i=0;i<G->n-1;i++) encoded[i]=cart_of_pol(G,pol_encoded[i]);
  free(pol_encoded);
  free(pol_msg);
}

void evaluate_rs_pol(u_int8_t* msg, int msg_len, u_int8_t* g, galois* G, u_int8_t* encoded) {
  int encoded_len = G->n-1;
  int deg_g = encoded_len-msg_len;
  u_int8_t* c_msg = (u_int8_t*)calloc(encoded_len, sizeof(u_int8_t));
  u_int8_t* r = (u_int8_t*)calloc(encoded_len, sizeof(u_int8_t));
  for(int i=deg_g; i<encoded_len; i++) {
    c_msg[i]=msg[i-deg_g];
    r[i] = c_msg[i];
  }

  //calcul du reste mod g
  for(int i=encoded_len-1; i>=deg_g; i--) {
    for(int j=i-deg_g; j<=i; j++) {
      u_int8_t m = (G->mult_table[r[i]*G->n + g[j-i+deg_g]]);
      // printf("i:%d j:%d r[j]:%d j-i+2:%d m:%d\n",i,j,r[j],j-i+2,m);
      r[j] = G->add_table[r[j]*G->n + m];
    }
  }

  //adding r to msg
  for(int i=0; i<deg_g; i++) c_msg[i] = r[i];

  for(int i=0; i<encoded_len;i++) {
    encoded[i] = c_msg[i];
  }
  // print_poly(c_msg,encoded_len-1);

  free(r);
  free(c_msg);

}

u_int8_t** encode_rs_8(u_int8_t** msg){
  galois* G = generate_galois_8();

  u_int8_t racines[2] = {2,3};
  u_int8_t* g = get_g(G, racines, 2);


  u_int8_t** send_msg = (u_int8_t**) malloc(7*sizeof(u_int8_t*));
  evaluate_rs_old(msg,g,G,send_msg);


  free(G);
  free(g);
  return send_msg;
}

// int main(){

//   u_int8_t msg[5][3] = {{1,0,0},{0,1,1},{1,0,1},{0,0,0},{1,1,0}};
//   printf("message : [");
//   for(int i=0;i<5;i++){
//     printf("[%d,%d,%d]",msg[i][0],msg[i][1],msg[i][2]);
//   }
//   printf("]\n");
//   printf("----> Reed Solomon\n");
//   u_int8_t** encoded = encode_rs_8(msg);
//   printf("encoded : [");
//   for(int i=0;i<7;i++){
//     printf("[%d,%d,%d]",encoded[i][0],encoded[i][1],encoded[i][2]);
//   }
//   printf("]\n");


//   return 0;

// }
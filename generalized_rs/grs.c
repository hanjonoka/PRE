#include "grs.h"

void evaluate_grs_old(u_int8_t** msg, u_int8_t V[7], u_int8_t lambda[7], galois* G, u_int8_t** encoded) {
  u_int8_t* f = (u_int8_t*)calloc(5,sizeof(u_int8_t));
  for (int i=0;i<5;i++){
    f[i]=pol_of_cart_8(G->elements,msg[i]);
  }

  u_int8_t* c_msg = (u_int8_t*)calloc(7,sizeof(u_int8_t));
  for(int i=0; i<7; i++) {
    // printf("i : %d",i);
    for(int j=0; j<5; j++) {
      // printf("\tj : %d",j);
      u_int8_t x_e_j = puiss_galois_8(G,V[i],j);
      // printf("\tV[i]^j : %d^%d = %d\n",V[i],j,x_e_j);

      u_int8_t monome = (G->mult_table[f[j]*8 + x_e_j]);
      // printf("\t\tf[j] * V[i]^j : %d * (%d^%d) = %d\n",f[j],V[i],j,monome);

      c_msg[i] = G->add_table[c_msg[i]*8 + monome];
      // printf("\t\tcmsg[i]:%d\n",c_msg[i]);
    }
    // printf("%d\n",c_msg[i]);
    c_msg[i] = G->mult_table[c_msg[i]*8 + lambda[i]];
    // printf("\n");
  }

  for(int i=0;i<7;i++){
    int c = c_msg[i];
    encoded[i]=cart_of_pol(G, c);
  }

  free(c_msg);
}

u_int8_t* get_lambda(galois* G, int i) {
  int len = G->n-1;
  u_int8_t* lambda = (u_int8_t*)calloc(len,sizeof(u_int8_t));
  for (int k=0;k<len;k++) {
    lambda[k] = puiss_galois(G,2,k*i);
  }
  return lambda;
}

void evaluate_grs(u_int8_t** msg, int msg_len, u_int8_t* V, u_int8_t* lambda, galois* G, u_int8_t** encoded) {
  u_int8_t* pol_msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t));
  for(int i=0;i<msg_len;i++) pol_msg[i]=pol_of_cart(G,msg[i]);
  u_int8_t* pol_encoded = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
  evaluate_grs_pol(pol_msg,msg_len,V,lambda,G,pol_encoded);
  for(int i=0;i<G->n-1;i++) encoded[i]=cart_of_pol(G,pol_encoded[i]);
  free(pol_encoded);
  free(pol_msg);
}

void evaluate_grs_pol(u_int8_t* msg, int msg_len, u_int8_t* V, u_int8_t* lambda, galois* G, u_int8_t* encoded) {
  int encoded_len = G->n-1;
  for(int i=0;i<encoded_len;i++){
    encoded[i]=0;
  }


  for(int i=0; i<encoded_len; i++) {
    for(int j=0; j<msg_len; j++) {
      u_int8_t x_e_j = puiss_galois_8(G,V[i],j);
      // printf("\tV[i]^j : %d^%d = %d\n",V[i],j,x_e_j);

      u_int8_t monome = (G->mult_table[msg[j]*G->n + x_e_j]);
      // printf("\t\tf[j] * V[i]^j : %d * (%d^%d) = %d\n",f[j],V[i],j,monome);

      encoded[i] = G->add_table[encoded[i]*G->n + monome];
    }
    encoded[i] = G->mult_table[encoded[i]*G->n + lambda[i]];
  }

}

u_int8_t** encode_grs_8(u_int8_t** msg) {
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda[7] = {1,1,1,1,1,1,1}; //coeficients de normalisation.
  galois* G = generate_galois_8();
  u_int8_t** send_msg = (u_int8_t**)malloc(7*sizeof(u_int8_t*));
  evaluate_grs_old(msg,V,lambda,G,send_msg);
  return send_msg;
}


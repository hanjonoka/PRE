#include "grs.h"

void evaluate_grs_old(u_int8_t** msg, u_int8_t V[7], u_int8_t lambda[7], gallois* G, u_int8_t** encoded) {
  u_int8_t* f = (u_int8_t*)calloc(5,sizeof(u_int8_t));
  for (int i=0;i<5;i++){
    f[i]=pol_of_cart_8(G->elements,msg[i]);
  }

  u_int8_t* c_msg = (u_int8_t*)calloc(7,sizeof(u_int8_t));
  for(int i=0; i<7; i++) {
    // printf("i : %d",i);
    for(int j=0; j<5; j++) {
      // printf("\tj : %d",j);
      u_int8_t x_e_j = puiss_gallois_8(G,V[i],j);
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
    encoded[i]=cart_of_pol(G->elements, c);
  }

  free(c_msg);
}

void evaluate_grs(u_int8_t** msg, int msg_len, u_int8_t* V, u_int8_t* lambda, gallois* G, u_int8_t** encoded) {
  int encoded_len = G->n-1;
  u_int8_t* f = (u_int8_t*) calloc(msg_len, sizeof(u_int8_t));
  for (int i=0;i<msg_len;i++) {
    f[i]=pol_of_cart(G,msg[i]);
  }

  u_int8_t* c_msg = (u_int8_t*) calloc(encoded_len,sizeof(u_int8_t));
  for(int i=0; i<encoded_len; i++) {
    for(int j=0; j<msg_len; j++) {
      u_int8_t x_e_j = puiss_gallois_8(G,V[i],j);
      // printf("\tV[i]^j : %d^%d = %d\n",V[i],j,x_e_j);

      u_int8_t monome = (G->mult_table[f[j]*G->n + x_e_j]);
      // printf("\t\tf[j] * V[i]^j : %d * (%d^%d) = %d\n",f[j],V[i],j,monome);

      c_msg[i] = G->add_table[c_msg[i]*G->n + monome];
    }
  }

  for(int i=0;i<encoded_len;i++){
    int c = c_msg[i];
    encoded[i]=cart_of_pol(G->elements, c);
  }

  free(c_msg);

}

u_int8_t** encode_grs_8(u_int8_t** msg) {
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t lambda[7] = {1,1,1,1,1,1,1}; //coeficients de normalisation.
  gallois* G = generate_gallois_8();
  u_int8_t** send_msg = (u_int8_t**)malloc(7*sizeof(u_int8_t*));
  evaluate_grs_old(msg,V,lambda,G,send_msg);
  return send_msg;
}


#include "../corps_gallois/gallois.h"

//tableau de formes polaires de G8 (u_int8_t)
u_int8_t* get_g(gallois* G, u_int8_t* elts, int deg){
  u_int8_t* g = calloc(deg+1,sizeof(u_int8_t));
  g[0]=1;
  u_int8_t* tmp = calloc(deg+1,sizeof(u_int8_t));
  for(int i=0;i<deg;i++){
    for(int j=0;j<=deg;j++){
      tmp[j] = G->add_table[G->mult_table[g[j]*8 + (elts[i])]*8 + (j>0 ? g[j-1] : 0)];
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

u_int8_t** encode_rs_8(u_int8_t msg[5][3]){
  gallois* G = generate_gallois_8();
  print_gallois_8(G);

  u_int8_t racines[2] = {2,3};
  u_int8_t* g = get_g(G, racines, 2);

  print_poly(g,2);

  u_int8_t* c_msg = (u_int8_t*)calloc(7,sizeof(u_int8_t));
  u_int8_t* r = (u_int8_t*)calloc(7,sizeof(u_int8_t));
  for(int i=2;i<7;i++) {
    c_msg[i]=pol_of_cart_8(G->elements,msg[i-2]);
    r[i]=c_msg[i];
  }

  //calcul du reste mod g
  for(int i=6; i>=2; i--){
    for(int j=i-2; j<=i; j++){
      u_int8_t m = (G->mult_table[r[i]*8 + g[j-i+2]]);
      // printf("i:%d j:%d r[j]:%d j-i+2:%d m:%d\n",i,j,r[j],j-i+2,m);
      r[j] = G->add_table[r[j]*8 + m];
    }
  }
  print_poly(c_msg,6);
  print_poly(r,6);

  c_msg[0]=r[0];
  c_msg[1]=r[1];
  print_poly(c_msg,6);

  u_int8_t** send_msg = malloc(7*sizeof(u_int8_t*));
  for(int i=0;i<7;i++){
    int c = c_msg[i];
    send_msg[i]=cart_of_pol(G->elements, c);
  }


  free(G);
  free(c_msg);
  free(g);
  free(r);
  return send_msg;
}

int main(){

  u_int8_t msg[5][3] = {{1,0,0},{0,1,1},{1,0,1},{0,0,0},{1,1,0}};
  u_int8_t** encoded = encode_rs_8(msg);
  printf("[");
  for(int i=0;i<7;i++){
    printf("[%d,%d,%d]",encoded[i][0],encoded[i][1],encoded[i][2]);
  }
  printf("]\n");


  return 0;

}
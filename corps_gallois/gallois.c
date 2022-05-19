#include "gallois.h"

void mult_gallois_8(u_int8_t* a, u_int8_t* b, u_int8_t* res) {
  u_int8_t* c = (u_int8_t*)(6,sizeof(u_int8_t));
  for(int i=0;i<3;i++) {
    for( int j=0;i<3;i++) {
      c[i+j] = c[i+j] & (a[i] ^ b[j]);
    }
  }

  //calcul de res = c mod P
  for(int i=0;i<3;i++){
    res[i] = c[i];
  }
  if(c[3]){ //ajoute x^3 mod P
    res[0] = res[0]^1;
    res[2] = res[2]^1;
  }
  if(c[4]){ //ajoute x^4 mod P
    res[0] = res[0]^1;
    res[1] = res[1]^1;
    res[2] = res[2]^1;
  }
  if(c[5]){ //ajoute x^4 mod P
    res[0] = res[0]^1;
    res[1] = res[1]^1;
  }
  free(c);
}

void add_gallois_8(u_int8_t* a, u_int8_t* b, u_int8_t* res) {
  for(int i=0;i<3;i++) {
    res[i]=a[i]^b[i];
  }
}

// u_int8_t opp_pol_8(u_int8_t i) {
//   return (8-i) % 8;
// }

gallois* generate_gallois_8() {
  u_int8_t elements[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {1,1,0}, {0,1,1}};
  gallois* G = (gallois*)malloc(sizeof(gallois));
  G->elements = (u_int8_t**)calloc(8,sizeof(u_int8_t*));
  for(int i=0;i<8;i++){
    G->elements[i] = (u_int8_t*)calloc(3,sizeof(u_int8_t));
    G->elements[i][0] = elements[i][0];
    G->elements[i][1] = elements[i][1];
    G->elements[i][2] = elements[i][2];
  }
  G->add_table = (u_int8_t*)calloc(64, sizeof(u_int8_t));
  G->mult_table = (u_int8_t*)calloc(64, sizeof(u_int8_t));

  u_int8_t* e = (u_int8_t*)calloc(3,sizeof(u_int8_t));
  for(int i=0;i<8;i++) {
    for(int j=0;j<8;j++){
      add_gallois_8(G->elements[i],G->elements[j],e);
      G->add_table[i*8 + j] = pol_of_cart_8(G->elements,e);

      G->mult_table[i*8 + j] = (i*j) % 8;
    }
  }
  free(e);
  return G;
}

u_int8_t* cart_of_pol(u_int8_t** elts, int i){
  u_int8_t* res = (u_int8_t*)calloc(3, sizeof(u_int8_t));
  for(int j=0;j<3;j++){
    res[j]=elts[i][j];
  }
  return res;
};

int el_eq(u_int8_t a[3], u_int8_t b[3]){
  return (a[0]==b[0])&&(a[1]==b[1])&&(a[2]==b[2]) ? 1 : 0;
}
u_int8_t pol_of_cart_8(u_int8_t** elts, u_int8_t* e){
  for(int i=0;i<8;i++){
    if(el_eq(elts[i],e)) return i;
  }
  return 9;
}

void print_gallois_8(gallois* G){
  printf("[");
  for(int i=0;i<8;i++){
    printf("[%d,%d,%d]",G->elements[i][0],G->elements[i][1],G->elements[i][2]);
  }
  printf("]\n");

  printf("+   N 0 1 2 3 4 5 6\t*   N 0 1 2 3 4 5 6\n");
  for(int i=0; i<8; i++){
    i ? printf("%d : ",i-1) : printf("N : ");
    for(int j=0; j<8;j++){
      G->add_table[i*8+j] ? printf("%d ",G->add_table[i*8+j]-1) : printf("N ");
    }
    i ? printf("\t%d : ",i-1) : printf("\tN : ");
    for(int j=0; j<8;j++){
      G->mult_table[i*8+j] ? printf("%d ",G->mult_table[i*8+j]-1) : printf("N ");
    }
    printf("\n");
  }
  printf("\n");
}


// int main(){
//   u_int8_t tb[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {1,1,0}, {0,1,1}};
//   gallois* G = generate_gallois_8(tb);
//   print_gallois_8(G);
//   return 0;
// }


#include "gallois.h"


void mult_galois_8(u_int8_t* a, u_int8_t* b, u_int8_t* res) { //jsp si ça marche
  u_int8_t* c = (u_int8_t*)calloc(6,sizeof(u_int8_t));
  for(int i=0;i<3;i++) {
    for( int j=0;i<3;i++) {
      c[i+j] = c[i+j] ^ (a[i] & b[j]);
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
  if(c[5]){ //ajoute x^5 mod P
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

void mult_mod(u_int8_t* a, u_int8_t* b, u_int8_t* res, u_int8_t* P, int deg_P) { //pas fini
  u_int8_t* c = (u_int8_t*) calloc(2*deg_P,sizeof(u_int8_t));
  for(int i=0;i<deg_P;i++) {
    for(int j=0;j<deg_P;j++) {
      c[i+j] = c[i+j] ^ (a[i] & b[j]); 
    }
  }

  for(int i=0;i<deg_P;i++) {
    res[i] = c[i];
  }
  for(int i=deg_P;i<2*deg_P;i++) {

  }
}

void mod_P(u_int8_t* res, u_int8_t* a, int deg_a, u_int8_t* P, int deg_P) {
  u_int8_t* b = (u_int8_t*) calloc(deg_a+1, sizeof(u_int8_t));
  for(int i=0;i<deg_a+1;i++) {
    b[i]=a[i];
  }

  for(int i=deg_a;i>=deg_P;i--) {
    for(int j=i-deg_P;j<=i;j++) {
      u_int8_t m = b[i] & P[j-i+deg_P];
      b[j] = b[j] ^ m;
    }
  }

  for(int i=0;i<deg_P;i++) {
    res[i] = b[i];
  }
  free(b);
}

void add(u_int8_t* a, u_int8_t* b, u_int8_t* res, int deg) {
  for(int i=0; i<deg+1; i++) {
    res[i]=a[i]^b[i];
  }
}

u_int8_t inverse(u_int8_t x, galois* G) {
  // printf("%d-(%d-1)%% %d\n",(G->n-1),x,(((G->n-1)-(x-1))%G->n-1)+1);
  return ((G->n-1)-(x-1)%(G->n-1))+1;
}

u_int8_t puiss_galois_8(galois* G, u_int8_t x, int n) {
  return n==0 ? 1 : G->mult_table[x*8 + puiss_galois_8(G, x, n-1)];
}

u_int8_t puiss_galois(galois* G, u_int8_t x, int n) {
  return n==0 ? 1 : G->mult_table[x*G->n + puiss_galois(G, x, n-1)];
}


galois* generate_galois_8() {
  u_int8_t elements[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {1,1,0}, {0,1,1}};
  galois* G = (galois*)malloc(sizeof(galois));
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

      G->mult_table[i*8 + j] = i==0 || j==0 ? 0 : ((i+j-2) % 7) + 1;
    }
  }
  free(e);
  return G;
}

galois* generate_galois(u_int8_t* P, int deg_P) {
  galois* G = (galois*)malloc(sizeof(galois));
  G->deg_P = deg_P;
  int n = 1<<(deg_P); //2^(deg_P+1) elements
  G->n = n;
  G->elements = (u_int8_t**)malloc(n * sizeof(u_int8_t*)); //2^(deg+1 elements)

  G->elements[0] = (u_int8_t*) calloc(deg_P, sizeof(u_int8_t));
  
  G->elements[1] = (u_int8_t*) calloc(deg_P, sizeof(u_int8_t));
  G->elements[1][0] = 1;

  u_int8_t* a = (u_int8_t*) calloc(n, sizeof(u_int8_t));
  for(int i=2; i<n; i++) {
    G->elements[i] = (u_int8_t*) calloc(deg_P, sizeof(u_int8_t));
    a[i-2] = 0;
    a[i-1] = 1;
    mod_P(G->elements[i], a, i-1, P, deg_P);
  }
  free(a);

  G->poids = (u_int8_t*) calloc(G->n,sizeof(u_int8_t));
  for(int i=0;i<G->n;i++) {
    for(int j=0;j<G->deg_P;j++) G->poids[i] += G->elements[i][j];
  }

  G->add_table = (u_int8_t*) calloc(n*n, sizeof(u_int8_t));
  G->mult_table = (u_int8_t*) calloc(n*n, sizeof(u_int8_t));

  u_int8_t* e = (u_int8_t*)calloc(deg_P,sizeof(u_int8_t));
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      add(G->elements[i], G->elements[j], e, deg_P-1);
      G->add_table[i*G->n + j] = pol_of_cart(G, e);

      G->mult_table[i*G->n + j] = i==0 || j==0 ? 0 : ((i+j-2) % (n-1)) + 1;
    }
  }
  free(e);

  return G;
}

u_int8_t* cart_of_pol(galois* G, int i){
  u_int8_t* res = (u_int8_t*)calloc(G->deg_P, sizeof(u_int8_t));
  for(int j=0;j<G->deg_P;j++){
    res[j]=G->elements[i][j];
  }
  return res;
};

int el_eq(u_int8_t* a, u_int8_t* b, int deg){
  for(int i=0; i<deg+1;i++) {
    if(a[i]!=b[i]) return 0;
  }
  return 1;
}
u_int8_t pol_of_cart_8(u_int8_t** elts, u_int8_t* e){
  for(int i=0;i<8;i++){
    if(el_eq(elts[i], e, 2)) return i;
  }
  printf("err pol_of_cart ");
  return -1;
}

u_int8_t pol_of_cart(galois* G, u_int8_t* e) {
  for(int i=0;i<G->n;i++) {
    if(el_eq(G->elements[i], e, G->deg_P-1)) return i;
  }
  return -1;
}

void print_galois_8(galois* G){
  printf("[");
  for(int i=0;i<8;i++){
    printf("[%d,%d,%d]",G->elements[i][0],G->elements[i][1],G->elements[i][2]);
  }
  printf("] nb=%d\n",G->n);

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

void print_galois(galois* G) {
  printf("[");
  for(int i=0;i<G->n;i++){
    printf("[%d",G->elements[i][0]);
    for(int j=1;j<G->deg_P;j++) {
      printf(",%d",G->elements[i][j]);
    }
    printf("]");
  }
  printf("] nb=%d\n",G->n);

  printf("[%d",G->poids[0]);
  for(int i=1;i<G->n;i++) {
    printf(",%d",G->poids[i]);
  }
  printf("]\n");

  printf(" +    N");
  for(int i=1;i<G->n;i++) i<=10 ? printf("  %d",i-1) : printf(" %d",i-1);
  printf("\t *    N");
  for(int i=1;i<G->n;i++) i<=10 ? printf("  %d",i-1) : printf(" %d",i-1);
  printf("\n\n");

  for(int i=0; i<G->n; i++){
    if(i==0) {printf(" N : ");}
    else if(i<=10) {printf(" %d : ",i-1);}
    else {printf("%d : ",i-1);}

    for(int j=0; j<G->n;j++){
      int k = G->add_table[i*G->n+j];
      if(k==0) {printf(" N ");}
      else if(k<=10) {printf(" %d ",k-1);}
      else {printf("%d ",k-1);}
    }

    printf("\t");
    if(i==0) {printf(" N : ");}
    else if(i<=10) {printf(" %d : ",i-1);}
    else {printf("%d : ",i-1);}

    for(int j=0; j<G->n;j++){
      int k = G->mult_table[i*G->n+j];
      if(k==0) {printf(" N ");}
      else if(k<=10) {printf(" %d ",k-1);}
      else {printf("%d ",k-1);}
    }
    printf("\n");
  }
  printf("\n");
}

#include "reed_solomon/rs.h"
#include "corps_gallois/gallois.h"
#include "math.h"
#include "stdio.h"

int comp_poly(u_int8_t* p1, u_int8_t* p2, int deg, galois* G) {
  u_int8_t a = inverse(p1[0],G);
  for(int i=0; i<=deg; i++) {
    // printf("%d %d ",p2[deg-i],(G->mult_table[p1[i]*G->n + a]));
    if(p2[deg-i] != (G->mult_table[p1[i]*G->n + a])) return 1;
  }
  return 0;
}

int main() {
  u_int8_t P32[6] = {1,1,1,0,1,1};
  u_int8_t P16[5] =  {1,1,0,0,1};
  u_int8_t P8[4] = {1,0,1,1};
  galois* G = generate_galois(P32, 5);
  int msg_len = 7;

  //init RS
  u_int8_t** g_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  int deg_g = G->n-1-msg_len;
  u_int8_t* racines = calloc(deg_g,sizeof(u_int8_t));
  for(int i=0;i<G->n-1;i++) {
    for(int j=0;j<deg_g;j++) racines[j] = ((i+j)%(G->n-1)+1);
    g_tab[i] = get_g(G,racines,deg_g);
    printf("%d : ",i);print_poly(g_tab[i],deg_g);
  }

  for(int i=0;i<G->n-1;i++) {
    int j=0;
    for(j=0;j<G->n-1;j++) {
      // printf("(%d %d) : ",i,(i+j)%(G->n-1));
      if(!comp_poly(g_tab[i],g_tab[(i+j)%(G->n-1)],deg_g, G)) break;
    }

    printf("g_%d",i);
    j<G->n-1 ? printf(" ~ g_%d\n",(i+j)%(G->n-1)) : printf("\n");
  }


  free(racines);
  for(int i=0;i<G->n-1;i++) free(g_tab[i]);
  free(g_tab);
}
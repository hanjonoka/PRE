#include "reed_solomon/rs.h"
#include "generalized_rs/grs.h"
#include "corps_gallois/gallois.h"
#include "math.h"
#include "stdio.h"

void print_matrix(u_int8_t* M, int heigth, int width) {
  for(int i=0;i<heigth;i++) {
    printf("| ");
    for(int j=0;j<width;j++) {
      printf("%d ",M[i*width + j]);
    }
    printf("|\n");
  }
}

u_int8_t* pivot_gauss(u_int8_t* M, int heigth, int width, galois* G) {
  u_int8_t* P = (u_int8_t*) calloc(width*heigth,sizeof(u_int8_t));
  for(int i=0;i<width*heigth;i++) P[i]=M[i]; //on recopie pour ne pas modifier M

  for(int j=0;j<width && j<heigth;j++) {
    // printf("-------------------\n");
    // print_matrix(P,heigth,width);
    //vérifier que P[j,j]!=0 et trouver un pivot non nul
    int p = P[j*width + j];
    for(int k=0;k<width;k++) {
      P[j*width + k]=G->mult_table[P[j*width + k]*G->n + inverse(p,G)]; //division de la ligne par le pivot
    }
    // printf("division\n");
    // print_matrix(P,heigth,width);
    // printf("soustraction\n");
    for(int i=0;i<heigth;i++) {
      int a = P[i*width + j];
      if(i!=j) {
       for(int k=0;k<width;k++) P[i*width + k]=G->add_table[P[i*width + k]*G->n + G->mult_table[a*G->n + P[j*width+k]]]; //P[i,k]-=P[j,k]P[i,j]
      }
    }
  }
  return P;
}

int main() {
  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);

  //init RS
  u_int8_t** g_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int i=0;i<G->n-1;i++) {
    u_int8_t racines[2] = {i+1,(i+1)%(G->n-1)+1};
    g_tab[i] = get_g(G,racines,2);
    printf("(%d,%d) : ",i,(i+1)%(G->n-1));
    print_poly(g_tab[i],2);
  }

  //init GRS
  u_int8_t V[7] = {1,2,3,4,5,6,7}; //tableau de formes polaires d'elements de G distincts où sera évalué le polynome.
  u_int8_t** lambda_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int i=0;i<G->n-1;i++) {
    lambda_tab[i] = get_lambda(G,i);
    printf("%d : ",i);
    print_lambda(lambda_tab[i], G->n-1);
  }

  u_int8_t** rs_matrix_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int k=0;k<G->n-1;k++) {
    rs_matrix_tab[k] = calloc((G->n-1) * 5, sizeof(u_int8_t));
    for(int i=0;i<5;i++) {
      for(int j=0;j<7;j++) {
        if(j<i) {
          rs_matrix_tab[k][i*(G->n-1) + j] = 0;
        }else if(j<i+3){
          rs_matrix_tab[k][i*(G->n-1) + j] = g_tab[k][j-i];
        }else{
          rs_matrix_tab[k][i*(G->n-1) + j] = 0;
        }
      }
    }
    printf("%d : \n",k);
    print_matrix(rs_matrix_tab[k],5,7);
    printf("\n");
    print_matrix(pivot_gauss(rs_matrix_tab[k],5,7,G),5,7);
    printf("-------------------------------\n");
  }

  u_int8_t** grs_matrix_tab = (u_int8_t**) calloc(G->n-1,sizeof(u_int8_t*));
  for(int k=0;k<G->n-1;k++) {
    grs_matrix_tab[k] = calloc((G->n-1) * 5, sizeof(u_int8_t));
    for(int i=0;i<5;i++) {
      for(int j=0;j<7;j++) {
        grs_matrix_tab[k][i*(G->n-1) + j] = G->mult_table[lambda_tab[k][j]*G->n + (puiss_galois(G,V[j],i))];
      }
    }
    printf("%d : \n",k);
    print_matrix(grs_matrix_tab[k],5,7);
    printf("\n");
    print_matrix(pivot_gauss(grs_matrix_tab[k],5,7,G),5,7);
    printf("-------------------------------\n");
  }
  // print_matrix(rs_matrix_tab[1],5,7);
  // printf("\n");
  // print_matrix(pivot_gauss(rs_matrix_tab[1],5,7,G),5,7);
}

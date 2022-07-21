#include "../generalized_rs/grs.h"
#include "../corps_gallois/gallois.h"
#include "../decoder/decoder.h"
#include "math.h"
#include "stdio.h"
#include "time.h"

//retourne un nombre entre 0 et 1
double uniform_distrib() {
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

//retourne un nombre suivant la loi N(sigma, mu) avec la transformée de Box-Muller
double normal_distrib(float sigma, float mu) {
  double v1=uniform_distrib();
  double v2=uniform_distrib();
  double u = cos(2*3.14*v2)*sqrt(-2.*log(v1));
  return u*sigma+mu;
}

//retourne 0 si les nombres mots sont égaux
int compare_words(u_int8_t* a, u_int8_t* b, int len){
  for(int i=0; i<len; i++) {
    if(a[i]!=b[i]) return 1;
  }
  return 0;
}

int main() {
  double sigma = 0.5; // caractéristique du canal
  int nb_test = 10;
  srand(time(NULL));

  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  int msg_len = 5; //k

  u_int8_t* V = get_lambda(G,1);
  u_int8_t* lambda = get_lambda(G,0);

  u_int8_t* msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t)); //mot à transmettre
  u_int8_t* encoded = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t)); //mot de code
  double* received = (double*) calloc(G->n-1,sizeof(double)); //mot après distortion par le canal

  for(u_int8_t j=0; j<pow(G->n,msg_len); j++) {
    int n=j;
    for(int k = 0; k<msg_len; k++) {
      msg[k] = n%G->n;
      n=n/G->n;
    }

    printf("msg : ");print_word(msg, msg_len);
    evaluate_grs_pol(msg, msg_len, V, lambda, G, encoded);
    printf("encoded : ");print_word(encoded, G->n-1);

    for(int k=0; k<nb_test; k++) { //on génère nb_test mots déformés par le canal
      for(int l=0; l<G->n-1; l++) {
        received[l] = normal_distrib(sigma, encoded[l]);
      }
      u_int8_t* decoded = decode_soft(received, G, msg_len, 4, sigma);
      printf("decoded : ");print_word(decoded, msg_len);
      if(compare_words(msg, decoded, msg_len)) {
        printf("erreur décodage\n");
      }else{
        printf("décodage réussi\n");
      }

      free(decoded);
    }

  }

}
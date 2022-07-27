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

  int nb_test = 1;

  srand(time(NULL));

  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  int msg_len = 5; //k
  int nb_msg = pow(G->n,msg_len);
  int taille_echantillon = nb_msg/50;

  u_int8_t* V = get_lambda(G,1);
  u_int8_t* lambda = get_lambda(G,0);


  FILE* courbe_taux_sigma = fopen("courbe_taux_sigma.csv","w");
  fprintf(courbe_taux_sigma,"sep=,\n");
  fprintf(courbe_taux_sigma,"sigma,max_mult,taux_reussite\n");
  fflush(courbe_taux_sigma);

  int nb_mult=1;
  for(int c=0; c<nb_mult; c++) {
    int max_mult = 6+c; //caractéristique du decoder

    int nb_sigma=8;
    float pas_sigma=0.0625;
    #pragma omp parallel for
    for(int i=1;i<nb_sigma;i++) {

      double sigma=(i+1)*pas_sigma; //caractéristique du canal
      int nb_reussite=0;
      u_int8_t* msg = (u_int8_t*) calloc(msg_len,sizeof(u_int8_t)); //mot à transmettre
      u_int8_t* encoded = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t)); //mot de code
      double* received = (double*) calloc(G->n-1,sizeof(double)); //mot après distortion par le canal

      for(int j=0; j<taille_echantillon; j++) {
        int n=rand()%nb_msg; //on selectionne des messages au hasard parmis l'ensemble pour avoir un échantillon représentatif (10%)
        for(int k = 0; k<msg_len; k++) {
          msg[k] = n%G->n;
          n=n/G->n;
        }

        // printf("msg : ");print_word(msg, msg_len);
        evaluate_grs_pol(msg, msg_len, V, lambda, G, encoded);
        // printf("encoded : ");print_word(encoded, G->n-1);

        for(int k=0; k<nb_test; k++) { //on génère nb_test mots déformés par le canal
          for(int l=0; l<G->n-1; l++) {
            received[l] = normal_distrib(sigma, encoded[l]);
          }
          u_int8_t* decoded = decode_soft(received, G, msg_len, max_mult, sigma);
          // printf("decoded : ");print_word(decoded, msg_len);
          if(!compare_words(msg, decoded, msg_len)) {
            nb_reussite++;
          }
          free(decoded);
        }
        if(j%10==0) {
          // print_progressbar((i*taille_echantillon)+(j+1),taille_echantillon*nb_sigma);
          printf("i=%d S = %.2f%% j=%d/%d max_mult=%d sigma=%.2f\n",i,((double)nb_reussite)/(((double)(j+1.))*((double)nb_test))*100.,(j+1),taille_echantillon,max_mult,sigma); fflush(stdout);
        }
      }
      double taux_reussite = ((double)nb_reussite)/((double)taille_echantillon)*100.;
      fprintf(courbe_taux_sigma,"%lf,%d,%lf\n",sigma,max_mult,taux_reussite);
      fflush(courbe_taux_sigma);
    }
  }
  fclose(courbe_taux_sigma);
}
#include "decoder.h"
#include "string.h"
#include "../generalized_rs/grs.h"

void print_reliability_matrix(double* M, int height, int width) {
  printf("PI =\n");
  for(int i=0; i<height; i++) {
    for(int j=0; j<width; j++) {
      printf("%lf ", M[i*width + j]);
    }
    printf("\n");
  }
}

void print_multiplicity_matrix(u_int8_t* M, int height, int width) {
  printf("M =\n");
  for(int i=0; i<height; i++) {
    for(int j=0; j<width; j++) {
      printf("%d ", M[i*width + j]);
    }
    printf("\n");
  }
}

void print_poly(u_int8_t** Q, int degx, int degy) {
  printf("Q = ");
  for(int i=0;i<=degx;i++) {
    for(int j=0; j<=degy; j++) {
      if(Q[i][j]!=0) printf(" + %d*x^%d*y^%d", Q[i][j],i,j);
    }
  }
  printf("\n");
}

void print_f(u_int8_t* f, int deg) {
  printf("f(x) = %d",f[0]);
  for(int i=1; i<=deg; i++) {
    printf(" + %d*x^%d", f[i],i);
  }
  printf("\n");
}

void print_list_f(u_int8_t** l_f, int l, int deg) {
  printf("list of f =\n");
  for(int i=0; i<l; i++) {
    printf("\t");
    print_f(l_f[i],deg);
  }
}

void print_poly_x_0(u_int8_t** Q, int degx, int degy) {
  printf("Q = ");
  for(int j=0; j<degy; j++) {
    if(Q[0][j]!=0) printf(" + %d*y^%d", Q[0][j],j);
  }
  printf("\n");
}

void print_word(u_int8_t* w, int l) {
  printf("w = [%d",w[0]);
  for(int i=1; i<l; i++) {
    printf(",%d",w[i]);
  }
  printf("]\n");
}

void print_list_word(u_int8_t** lw, int n, int l) {
  printf("List of words = \n");
  for(int i=0; i<n; i++) {
    printf("\t");
    print_word(lw[i],l);
  }
}

void print_progressbar(int progress, int total) {
  printf("\rProgress : ");
  for(int i=0; i<progress; i++) {
    printf("#");
  }
  for(int i=0; i<total-progress; i++) {
    printf("-");
  }
  printf(" %d%%",progress*100/total);
  fflush(stdout);
}

//renvoie un corps de galois de polynôme générateur de degré d.
galois* init_galois(int d) {
  if(d==1) {
    u_int8_t P[2] = {1,1};
    galois* G = generate_galois(P, 1);
    return G;
  }else if(d==2) {
    u_int8_t P[3] =  {1,1,1};
    galois* G = generate_galois(P, 2);
    return G;
  }else if(d==3) {
    u_int8_t P[4] =  {1,0,1,1};
    galois* G = generate_galois(P, 3);
    return G;
  }else if(d==4) {
    u_int8_t P[5] =  {1,1,0,0,1};
    galois* G = generate_galois(P, 4);
    return G;
  }else if(d==5) {
    u_int8_t P[6] = {1,1,1,0,1,1};
    galois* G = generate_galois(P, 5);
    return G;
  }else{
    printf("Impossible de créer un corps de galois avec d=%d\n",d);
    return NULL;
  }
}


int main(int argc, char* argv[]) {
  galois* G;
  int k;
  u_int8_t* msg;
  u_int8_t* V;
  u_int8_t* lambda;
  double* PI;
  int i_opt;
  int max_mult=4;

  if (argc==1){
    printf("Pas assez d'arguments. --help pour l'aide\n");
    return 1;
  }else if(argc>1 && !strcmp(argv[1], "--help")) {
    printf("Utilisation : %s d k (-c mot hard|soft)|(-f mot)\n", argv[0]);
    printf("\td : degré du polynôme générateur du corps de Galois (1<=d<=6)\n");
    printf("\tk : nombre de symbole dans le mot à encoder\n");
    printf("\tmot sous la forme : x1,...,xm\n");
    printf("Options :\n");
    printf("-c : mot à décoder. Dans ce cas m=n\n");
    printf("-f : mot à encoder puis décoder. Dans ce cas m=k\n");
    printf("-max_mult : multiplicité maximale dans la matrice\n");
    return 0;
  }else if(argc>=5) {
    int d = atoi(argv[1]);
    G = init_galois(d);
    if(G==NULL) return 1;
    V = get_lambda(G,1);
    lambda = get_lambda(G,0);


    k=atoi(argv[2]);
    if(k<=1 || k>G->n-1) {
      printf("k doit être entre 2 et n\n");
      return 0;
    }

    if(!strcmp(argv[3],"-c")) {
      if(!strcmp(argv[5],"hard")) {
        msg = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
        char* ptr = strtok(argv[4],",");
        int i=0;
        while(ptr!=NULL) {
          msg[i]=atoi(ptr);
          ptr = strtok(NULL,",");
          i++;
        }
        if(i!=G->n-1) {
          printf("Le mot doit comporter n caratères\n");
          return 1;
        }
        PI = generate_reliabitity_matrix_hard(G, msg);
      }else if(!strcmp(argv[5],"soft")) {
        printf("not yet implemented.\n");
        return 1;
      }else{
        printf("%s : Soft ou Hard ??\n",argv[5]);
        return 1;
      }
      i_opt = 6;
    } else if(!strcmp(argv[3],"-f")) {
      u_int8_t* f = (u_int8_t*) calloc(k,sizeof(u_int8_t));
      msg = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));

      char* ptr = strtok(argv[4],",");
      int i=0;
      while(ptr!=NULL) {
        f[i]=atoi(ptr);
        ptr = strtok(NULL,",");
        i++;
      }
      if(i!=k) {
        printf("Le mot doit comporter k caratères\n");
        return 1;
      }

      printf("Mot à encoder : ");
      print_word(f,k);
      evaluate_grs_pol(f,k,V,lambda,G,msg);
      PI = generate_reliabitity_matrix_hard(G, msg);

      i_opt=5;
    }

    for(int i=i_opt; i<argc; i++) {
      if(!strncmp(argv[i],"-max_mult=",10)) {
        max_mult=atoi(argv[i]+(10*sizeof(char)));
      }
    }
  }else{
    printf("Impossible d'interpréter les arguments. --help pour l'aide\n");
  }

  printf("Mot à décoder : ");
  print_word(msg,G->n-1);

  printf("\nMultiplicity Assignement\n");
  print_reliability_matrix(PI, G->n-1, G->n);
  u_int8_t* M = generate_multiplicity_proportional(G, PI, max_mult, k);
  printf("Result :\n");
  print_multiplicity_matrix(M,G->n-1,G->n);

  printf("\nInterpolation\n");
  u_int8_t** Q = interpolate(G, M, G->n-1, G->n, k);
  int omega = compute_omega(M, G->n-1, G->n, k);
  int L = omega/(k-1);
  int c = cost(M, G->n-1, G->n);
  printf("Result :\n");
  print_poly(Q,c,L);

  printf("\nFactorization\n");
  int n_w;
  int degx_norm, degy_norm;
  u_int8_t** Q_norm = normalize(G, Q, c, L, &degx_norm, &degy_norm);
  // print_poly(Q_norm, degx_norm, degy_norm);
  u_int8_t** l_w = Factorize(G, Q_norm, degx_norm, degy_norm, 0, k, &n_w);

  printf("Results :\n");
  print_list_f(l_w, n_w, k-1);
  u_int8_t** l_sent = (u_int8_t**) calloc(n_w,sizeof(u_int8_t*));
  if(n_w==0) {
    printf("Le décodage à échoué.\n");
    return 1;
  }
  for(int i=0;i<n_w;i++) {
    l_sent[i] = (u_int8_t*) calloc(G->n-1,sizeof(u_int8_t));
    evaluate_grs_pol(l_w[i],k,V,lambda,G,l_sent[i]);
    print_word(l_sent[i],G->n-1);
    printf("score = %d\n",score(l_sent[i],G->n-1, M, G->n-1, G->n));
  }


  printf("\nSelection du meilleur mot\n");
  int score_max = score(l_sent[0],G->n-1, M, G->n-1, G->n);
  int i_max = 0;
  for(int i=1; i<n_w; i++) {
    int s = score(l_sent[i],G->n-1, M, G->n-1, G->n);
    if(score_max<s) {
      score_max=s;
      i_max=i;
    }
  }
  print_word(l_sent[i_max],G->n-1);
  printf("score = %d\n",score_max);
  print_f(l_w[i_max], k-1);
  return 0;
}
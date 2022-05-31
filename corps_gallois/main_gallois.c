#include "gallois.h"

int main(){
  u_int8_t P[4] =  {1,0,1,1};
  galois* G = generate_galois(P, 3);
  print_galois(G);

  u_int8_t P2[5] =  {1,1,0,0,1};
  galois* G2 = generate_galois(P2, 4);
  print_galois(G2);
  return 0;
}

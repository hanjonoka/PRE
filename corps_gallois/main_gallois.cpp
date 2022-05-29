#include "gallois.h"

int main(){
  u_int8_t P[4] =  {1,0,1,1};
  gallois* G = generate_gallois(P, 3);
  print_gallois(G);

  u_int8_t P2[5] =  {1,1,0,0,1};
  gallois* G2 = generate_gallois(P2, 4);
  print_gallois(G2);
  return 0;
}

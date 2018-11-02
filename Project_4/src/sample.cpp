#include "sample.h"
// #include <chrono>
#include <time.h>
#include <cmath>
#include <iostream>

using namespace std;

sample::sample(int cycs, int d, double T){
  temperature = T;
  dim = 2;
  meanE = meanE2 = meanM = meanM2 = absM = 0;
  nCycles = cycs;

}

void sample::initializeLattice(int ordered){
  srand(clock());

  //allocate memory
  lattice = new double*[dim+2];
  for(int i = 0; i < dim+2; i++) lattice[i] = new double[dim+2];

  //set all lattice values to "spin up"
  if (ordered == 1) {
    for(int i = 0; i < dim+2; i++){
      for(int j = 0; j < dim+2; j++){
        lattice[i][j] = 1;
      }
    }
    meanM = dim*dim;
  }//end ordered lattice

  //set random lattice values
  else {
    for(int i = 0; i < dim+2; i++){
      for(int j = 0; j < dim+2; j++){
        double randNum = (((double) rand() / (RAND_MAX)) > 0.5) ? 1: -1;
        lattice[i][j] = randNum;
        meanM += (double) lattice[i][j];
      }
    }

    //set ghost vector values
    for(int i = 0; i < dim+1; i++){
      lattice[i][dim+1] = lattice[i][1];
      lattice[i][0] = lattice[i][dim];
      lattice[dim+1][i] = lattice[1][i];
      lattice[0][i] = lattice[dim][i];
    }
  }//end random lattice

  lattice[0][0] = lattice[0][dim+1] = lattice[dim+1][0] = lattice[dim+1][dim+1] = 0;

  //initial energy
  for(int i = 1; i < dim+1; i++){
    for(int j = 1; j < dim+1; j++){
      meanE -= lattice[i][j]*(lattice[i+1][j] + lattice[i][j+1]);
    }
  }
}//end initialize lattice function


void sample::precalculateProbs(){
  for(int i = 0; i < 17; i+=4) this->w[i] = exp(-(i-8)/(this->temperature));
}


void sample::metropolis(){
  srand(clock());
  int x, y, r;
  double E, M = 0;

  for(int cycle = 0; cycle < nCycles; cycle++){

    for(int i = 0; i < dim*dim+1; i++){
      //pick random element in lattice
      x = 1 + (rand() & (int)(dim-1));
      y = 1 + (rand() & (int)(dim-1));

      int dE = 2*lattice[x][y]*(lattice[x-1][y] + lattice[x+1][y] + //calculate delta E
        lattice[x][y-1] + lattice[x][y+1]);
        r = (double) rand() / (RAND_MAX);
        if(r <= w[dE+8]){
          lattice[x][y] *= -1;
          if(x == 1) lattice[dim+1][y] *= -1;
          if(x == dim) lattice[0][y] *= -1;
          if(y == 1) lattice[x][dim+1] *= -1;
          if(y == dim) lattice[x][0] *= -1;
          E += dE;
          M += 2*lattice[x][y];
        }
    }//end lattice sweep

    meanE += E;
    meanE2 += E*E;
    meanM += M;
    meanM2 += M*M;
    absM += fabs(M);
  }//end MC loop
}//end metropolis function


void sample::normalize(){
  meanE /= dim*dim*nCycles;
  meanE2 /= dim*dim*nCycles;
  meanM /= dim*dim*nCycles;
  meanM2 /= dim*dim*nCycles;
  absM /= dim*dim*nCycles;
}

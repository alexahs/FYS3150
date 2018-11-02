#include "sample.h"
// #include <chrono>
#include <time.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

sample::sample(int cycs, int d, double T, int start, int stop){
  temperature = T;
  dim = d;
  nCycles = cycs;
  meanE = 0;
  meanE2 = 0;
  meanM = 0;
  meanM2 = 0;
  absM = 0;
  loopStart = start;
  loopStop = stop;


}

void sample::initializeLattice(int ordered){
  // srand(clock());
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  //allocate memory
  lattice = new double*[dim+2];
  for(int i = 0; i < dim+2; i++) lattice[i] = new double[dim+2];

  //set all lattice values to "spin up"
  if (ordered == 1) {
    for(int i = 1; i < dim+1; i++){
      for(int j = 1; j < dim+1; j++){
        lattice[i][j] = 1;
      }
    }
    meanM = dim*dim;
  }//end ordered lattice

  //set random lattice values
  else {
    for(int i = 1; i < dim + 1; i++){
      for(int j = 1; j < dim + 1; j++){
        double randNum = (distribution(gen) > 0.5) ? 1: -1;
        lattice[i][j] = randNum;
        meanM += (double) lattice[i][j];
      }
    }
  }//end random lattice


  //set ghost vector values
  for(int i = 0; i < dim+2; i++){
    lattice[i][dim+1] = lattice[i][1];
    lattice[i][0] = lattice[i][dim];
    lattice[dim+1][i] = lattice[1][i];
    lattice[0][i] = lattice[dim][i];
  }

  lattice[0][0] = lattice[0][dim+1] = lattice[dim+1][0] = lattice[dim+1][dim+1] = 0;

  //initial energy
  for(int i = 1; i < dim+1; i++){
    for(int j = 1; j < dim+1; j++){
      meanE -= lattice[i][j]*(lattice[i+1][j] + lattice[i][j+1]);
    }
  }
}//end initialize lattice function


void sample::precalculateProbs(){
  for(int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
}


void sample::metropolis(){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // srand(clock());
  int x, y;
  double r;
  double E = 0;
  double M = 0;

  for(int cycle = loopStart; cycle < loopStop; cycle++){

    for(int i = 1; i < dim*dim+1; i++){
      //pick random element in lattice
      x = 1 + (int) (distribution(gen)*(double)dim);
      y = 1 + (int) (distribution(gen)*(double)dim);

      int dE = 2*lattice[x][y]*(lattice[x-1][y] + lattice[x+1][y] + //calculate delta E
        lattice[x][y-1] + lattice[x][y+1]);
        r = distribution(gen);
        if(r <= w[dE+8]){
          lattice[x][y] *= -1;
          if(x == 1) lattice[dim+1][y] *= -1;
          if(x == dim) lattice[0][y] *= -1;
          if(y == 1) lattice[x][dim+1] *= -1;
          if(y == dim) lattice[x][0] *= -1;
          E += (double) dE;
          M += 2*lattice[x][y];
        }
    }//end lattice sweep

    // if(cycle > 0.2*nCycles)
    meanE += E;
    meanE2 += E*E;
    meanM += M;
    meanM2 += M*M;
    absM += fabs(M);
  }//end MC loop
}//end metropolis function


void sample::normalize(){
  meanE /= dim*dim*nCycles + 1;
  meanE2 /= dim*dim*nCycles;
  meanM /= dim*dim*nCycles + 1;
  meanM2 /= dim*dim*nCycles;
  absM /= dim*dim*nCycles;
}

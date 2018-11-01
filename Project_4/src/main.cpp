#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <vector>
#include <string>
#include <time.h>
#include <stdlib.h>
#include "initializeLattice.h"
#include "analytical.h"


using namespace std;

void metropolis(int nCycles, int dim, double **lattice, double *expVals, double *w);
void printLattice(double **lattice, int dim);

int main(int argc, char* argv[])
{
  srand(clock());

  //problem parameters
  int dim = 2;
  double T = 2.0;

  int N = dim*dim;
  int nCycles = 1000000;
  double norm = 1.0/((double) nCycles * (double) N);

  //allocate memoey
  double *expectationValues = new double[5];
  double *w = new double[17];
  double **lattice = new double*[dim+2];
  for(int i = 0; i < dim+2; i++) lattice[i] = new double[dim+2];



  //initialize arrays values
  for(int i = 0; i < 5; i++) expectationValues[i] = 0.0;
  for(int i = 0; i < 17; i+=4) w[i] = exp(-(i-8)/T);


  randomLattice(lattice, dim, expectationValues);
  // orderedLattice(lattice, dim, expectationValues);

  printLattice(lattice, dim);



  cout << "Initial E: " << expectationValues[0]/((double) N) << endl;
  cout << "Initial E2: " << expectationValues[1]/((double) N) << endl;
  cout << "Initial M: " << expectationValues[2]/((double) N) << endl;
  cout << "Initial M2: " << expectationValues[3]/((double) N) << endl;
  cout << "Initial |M|: " << expectationValues[4]/((double) N) << endl;
  cout << "=======================" << endl;

  metropolis(nCycles, dim, lattice, expectationValues, w);

  cout << "Final E: " << expectationValues[0]*norm << endl;
  cout << "Final E2: " << expectationValues[1]*norm << endl;
  cout << "Final M: " << expectationValues[2]*norm << endl;
  cout << "Final M2: " << expectationValues[3]*norm << endl;
  cout << "Final |M|: " << expectationValues[4]*norm << endl;
  // cout << "Expected energy: " << energyExpected() << endl;
  // cout << "Expected magnetization: " << magnSquaredExpected() << endl;


  return 0;
}

void metropolis(int nCycles, int dim, double **lattice, double T0, double T1, int resolution){



  int nSpins = dim*dim;
  int x, y;
  double r;

  dT = T1 -T0;
  Tstep = 1.0/resolution;



  double *expectationValues = new double[5];
  double **lattice = new double*[dim+2];
  double *w = new double[17];
  for(int i = 0; i < 5; i++) expectationValues[i] = 0.0;
  for(int i = 0; i < dim+2; i++) lattice[i] = new double[dim+2];


  //temperature loop
  for(int t = 0; t < resolution; t++){


    T = T0 + dT*Tstep;
    
    for(int i = 0; i < 5; i++) expectationValues[i] = 0.0;
    randomLattice(lattice, dim, expectationValues);
    // orderedLattice(lattice, dim, expectationValues);


    for(int i = 0; i < 17; i+=4) w[i] = exp(-(i-8)/T);

    //metropolis
    for(int cycle = 0; cycle < nCycles; cycle++){
      //Flip spins
      for(int i = 1; i < nSpins+1; i++){
        //pick random element in lattice
        x = 1 + (rand() & (int)(dim-1));
        y = 1 + (rand() & (int)(dim-1));


        int dE = 2*lattice[x][y]*(lattice[x-1][y] + lattice[x+1][y] + //calculate delta E
          lattice[x][y-1] + lattice[x][y+1]);
          r = (double) rand() / (RAND_MAX);
          if(r <= w[dE+8]){ //Pick condition
            lattice[x][y] *= -1;
            if(x == 1) lattice[dim+1][y] *= -1;
            if(x == dim) lattice[0][y] *= -1;
            if(y == 1) lattice[x][dim+1] *= -1;
            if(y == dim) lattice[x][0] *= -1;
            E += dE;
            M += 2*lattice[x][y];

          }
        } //end lattice sweep
        expVals[0] += E;
        expVals[1] += E*E;
        expVals[2] += M;
        expVals[3] += M*M;
        expVals[4] += fabs(M);

    } //end MC loop
  } //end temperature loop
} //end metropolis function

void printLattice(double **lattice, int dim){
  //Prints the lattice to terminal, including boundary vectors
  cout << "Matrix:" << endl;
  for(int i = 0; i < dim+2; i++){
    for(int j = 0; j < dim+2; j++){
      cout << setw(12) << setprecision (10) << lattice[i][j];
    }
    cout << endl;
  }
}

void output(double T, double nCycles, double dim, double *expVals){
  cout << setw(12) << setprecision (10) << T;
  cout << setw(12) << setprecision (10) << expVals[0];
  cout << setw(12) << setprecision (10) << expVals[1];
  cout << setw(12) << setprecision (10) << expVals[2];
  cout << setw(12) << setprecision (10) << expVals[3];
  cout << setw(12) << setprecision (10) << expVals[4] << endl;
}

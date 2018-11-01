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

void metropolis(int nCycles, int dim, double **lattice, double *expVals, vector<double> w);
void printLattice(double **lattice, int dim);

int main(int argc, char* argv[])
{
  srand(clock());



  //problem parameters
  int dim = 2;
  double T = 1;
  int N = dim*dim;
  int nCycles = 1000000;
  double E, M, E2, M2 = 0;
  double *expectationValues = new double[5];
  // vector<double> expectationValues(5);

  //allocate memoey
  double **lattice = new double*[dim+2];
  vector<double> w(17);

  //initialize dEnergy values
  for(int i = 0; i < 5; i++) expectationValues[i] = 0.0;
  for(int i = 0; i < dim+2; i++) lattice[i] = new double[dim+2];
  for(int i = 0; i < 17; i+=4) w[i] = exp(-(i-8)/T);


  randomLattice(lattice, dim, expectationValues);
  // orderedLattice(lattice, dim, expectationValues);
  cout << "Initial E: " << expectationValues[0] << endl;
  cout << "Initial E2: " << expectationValues[1] << endl;
  cout << "Initial M: " << expectationValues[2] << endl;
  cout << "Initial M: " << expectationValues[3] << endl;
  cout << "Initial |M|: " << expectationValues[4] << endl;

  metropolis(nCycles, dim, lattice, expectationValues, w);


  //print initial matrix
  // printLattice(lattice, dim);

    //print final lattice
    // printLattice(lattice, dim);

  cout << "Final E: " << expectationValues[0] << endl;
  cout << "Final E2: " << expectationValues[1] << endl;
  cout << "Final M: " << expectationValues[2] << endl;
  cout << "Final M: " << expectationValues[3] << endl;
  cout << "Final |M|: " << expectationValues[4] << endl;
  // cout << "Expected energy: " << energyExpected() << endl;
  // cout << "Expected magnetization: " << magnExpected() << endl;


  return 0;
}

void metropolis(int nCycles, int dim, double **lattice, double *expVals, vector<double> w){

  int nSpins = dim*dim;
  int x, y;
  double r;
  double E = 0;
  double M = 0;

  for(int cycle = 0; cycle < nCycles; cycle++){
    //Flip spins
    for(int i = 1; i < nSpins+1; i++){
      //pick random element in lattice
      x = 1 + (rand() & (int)(dim-1));
      y = 1 + (rand() & (int)(dim-1));

      //calculate delta E
      int dE = 2*lattice[x][y]*(lattice[x-1][y] + lattice[x+1][y] +
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
      }
  }
  expVals[0] += E;
  expVals[1] += E*E;
  expVals[2] += M;
  expVals[3] += M*M;
  expVals[4] += fabs(M);
}

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

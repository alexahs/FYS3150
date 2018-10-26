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


using namespace std;

void randomLattice(double **lattice, int dim, double *E, double *M);
void orderedLattice(double **lattice, int dim, double *E, double *M);

int main(int argc, char* argv[])
{
  srand(clock());
  //allocate memoey
  int dim = 3;
  double T = 300;
  int N = dim*dim;
  int x, y;
  double E, M = 0;
  double r;
  double **lattice = new double*[dim+2];
  double w[17];
  for(int i = 0; i < dim+2; i++){
    lattice[i] = new double[dim+2];
  }
  for(int i = 0; i < 17; i+=4){
    w[i] = exp(-(i-8)/T);
  }

  //create lattice with random states
  randomLattice(lattice, dim, &E, &M);
  // orderedLattice(lattice, dim, &E, &M);
  cout << "Initial E: " << E << endl;
  cout << "Initial M: " << M << endl;

  //print initial matrix
  for(int i = 0; i < dim+2; i++){
    for(int j = 0; j < dim+2; j++){
      cout << setw(12) << setprecision (10) << lattice[i][j];
    }
    cout << endl;
  }
  cout << endl;


  //Flip states
  for(int i = 1; i < dim*dim+1; i++){
    //pick random element in lattice
    x = 1 + (rand() & (int)(dim-1));
    y = 1 + (rand() & (int)(dim-1));

    //calculate delta E
    int dE = 2*lattice[x][y]*(lattice[x-1][y] + lattice[x+1][y] +
      lattice[x][y-1] + lattice[x][y+1]);
      r = (double) rand() / (RAND_MAX);
      // cout << dE << endl;
      if(r <= w[dE+8]){ //Pick condition
        lattice[x][y] *= -1;
        if(x == 1) lattice[dim+1][y] *= -1;
        if(x == dim) lattice[0][y] *= -1;
        if(y == 1) lattice[x][dim+1] *= -1;
        if(y == dim) lattice[x][0] *= -1;
        E += dE;
        M += lattice[x][y]*2;
      }
    }

    for(int i = 0; i < dim+2; i++){
      for(int j = 0; j < dim+2; j++){
        cout << setw(12) << setprecision (10) << lattice[i][j];
      }
      cout << endl;
    }


  return 0;
}


void orderedLattice(double **lattice, int dim, double *E, double *M){
  for(int i = 0; i < dim+2; i++){
    for(int j = 0; j < dim+2; j++){
      lattice[i][j] = 1;
    }
  }
  lattice[0][0] = lattice[0][dim+1] = lattice[dim+1][0] = lattice[dim+1][dim+1] = 0;

  for(int i = 1; i < dim+1; i++){
    for(int j = 1; j < dim+1; j++){
      *M += lattice[i][j];
      *E -= lattice[i][j]*lattice[i+1][j] + lattice[i][j]*lattice[i][j+1];
    }
  }
}




void randomLattice(double **lattice, int dim, double *E, double *M){
  for(int i = 1; i < dim+1; i++){
    for(int j = 1; j < dim+1; j++){
      double randNum = (((double) rand() / (RAND_MAX)) > 0.5) ? 1: -1;
      lattice[i][j] = randNum;
    }
  }
  for(int i = 0; i < dim+1; i++){
    lattice[i][dim+1] = lattice[i][1];
    lattice[i][0] = lattice[i][dim];
    lattice[dim+1][i] = lattice[1][i];
    lattice[0][i] = lattice[dim][i];
  }

  for(int i = 1; i < dim+1; i++){
    for(int j = 1; j < dim+1; j++){
      *M += lattice[i][j];
      *E += lattice[i][j]*lattice[i+1][j] + lattice[i][j]*lattice[i][j+1];
    }
  }
}

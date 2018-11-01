#include "initializeLattice.h"
using namespace std;


void orderedLattice(double **lattice, int dim, double *expVals){
  //Initializes the lattice matrix with each element "spin up"
  double E_initial = 0;
  double M_initial = 0;

  for(int i = 0; i < dim+2; i++){
    for(int j = 0; j < dim+2; j++){
      lattice[i][j] = 1;
    }
  }
  lattice[0][0] = lattice[0][dim+1] = lattice[dim+1][0] = lattice[dim+1][dim+1] = 0;


  for(int i = 1; i < dim+1; i++){
    for(int j = 1; j < dim+1; j++){
      E_initial -= lattice[i][j]*lattice[i+1][j] + lattice[i][j]*lattice[i][j+1]; //E
      M_initial += lattice[i][j]; //M
    }
  }
  expVals[0] = E_initial;
  expVals[1] = E_initial*E_initial;
  expVals[2] = M_initial;
  expVals[3] = M_initial*M_initial;
  expVals[4] = fabs(M_initial);
}


void randomLattice(double **lattice, int dim, double *expVals){
  //Initializes the lattice matrix with random spin values
  double E_initial = 0;
  double M_initial = 0;

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
      E_initial += lattice[i][j]*lattice[i+1][j] + lattice[i][j]*lattice[i][j+1]; //E
      M_initial += lattice[i][j]; //M
    }
  }
  expVals[0] = E_initial;
  expVals[1] = E_initial*E_initial;
  expVals[2] = M_initial;
  expVals[3] = M_initial*M_initial;
  expVals[4] = fabs(M_initial);
}

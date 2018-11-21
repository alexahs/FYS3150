using namespace std;
#include <random>
#include <iostream>
#include "initializelattice.h"

void InitializeLattice ( int dim, double **SpinMatrix, double &E, double &M, int ordered) {
  //This function initializes the spin lattice. Arguments include the dimension length,
  //matrix of pointers, energy variable E and magnetization M, and an int 1/0 (ordered/random)
  // for specifying initial lattice configuration.

  random_device rd;
  mt19937_64 gen( rd() );
   // Uniform "REAL" distribution for x \in [0,1]
  uniform_real_distribution<double>RandomNumberGenerator( 0, 1 );


  //Sets up an ordered lattice with all spins up
  if(ordered == 1){
    for ( int x = 0; x < dim; x++ ) {
        for ( int y = 0; y < dim; y++ ) {
            SpinMatrix[x][y] = 1.0;
            M += ( double ) SpinMatrix[x][y];
        }
      }
    }

  else{
    //Sets up a lattice with random spins
    for ( int x = 0; x < dim; x++ ) {
        for ( int y = 0; y < dim; y++ ) {
          double randNum = (RandomNumberGenerator(gen) > 0.5) ? 1: -1;
          SpinMatrix[x][y] = randNum;
          M += ( double ) SpinMatrix[x][y];
        }
      }
    }
    //Initial energy of lattice
    for ( int x = 0; x < dim; x++ ) {
        for ( int y = 0; y < dim; y++ ) {
            E -= ( double ) SpinMatrix[y][x]*
                    ( SpinMatrix[PeriodicBoundary( y, dim, -1 )][x] +
                    SpinMatrix[y][PeriodicBoundary( x, dim, -1 )] );
        }
    }
}

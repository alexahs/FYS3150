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
#include <iomanip>
#include <mpi.h>
// #include "initializeLattice.h"
// #include "analytical.h"
#include "sample.h"
#include "analytical.h"
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using namespace std;


void printEnergies(sample Sample, double T, int dim, int nCycles, double *expVals);
void printMatrix(sample Sample, int dim);

int main(int argc, char *argv[]){
  int dim = 2;
  int ordered = 0;
  double temp_init = 1.0;
  double temp_final = 1.1;
  double temp_step = 0.01;
  int cycles = 1<<20;
  double norm = 1.0/(dim*dim);
  steady_clock::time_point programStart;
  programStart = steady_clock::now();

  int nProcs, my_rank;

  double *ExpecVals = new double[5];
  double *AllExpecVals = new double[5];
  for(int i = 0; i < 5; i++){
    ExpecVals[i] = 0;
    AllExpecVals[i] = 0;
  }


  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  MPI_Bcast (&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int cycleInterval = cycles/nProcs;
  int loopStart = my_rank*cycleInterval;
  int loopStop = (my_rank+1)*cycleInterval;

  for(double temperature = temp_init; temperature <= temp_final; temperature+=temp_step){

    sample S = sample(cycleInterval, dim, temperature, loopStart, loopStop);
    S.initializeLattice(ordered);
    S.precalculateProbs();
    S.metropolis();
    ExpecVals[0] = S.meanE;
    ExpecVals[1] = S.meanE2;
    ExpecVals[2] = S.meanM;
    ExpecVals[3] = S.meanM2;
    ExpecVals[4] = S.absM;
    for(int i = 0; i < 5; i++){
      MPI_Reduce(&ExpecVals[i], &AllExpecVals[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if(my_rank == 0){
      printEnergies(S, temperature, dim, cycles, AllExpecVals);
    }
  }

  if(my_rank == 0){
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);
    cout << "Program time: " << programTime.count() << endl;
  }


  MPI_Finalize();
  return 0;
}

void printEnergies(sample Sample, double T, int dim, int nCycles, double *expVals){
  for(int i = 0; i < 5; i++) expVals[i] /= (double) nCycles;
  double norm = 1.0/(dim*dim);
  cout << setw(15) << setprecision (6) << T;
  for(int i = 0; i < 5; i++){
    cout << setw(15) << setprecision (6) << expVals[i]*norm;
  }
  cout << endl;


  //
  // cout << setw(2) << setprecision (6) << Sample.temperature;
  // cout << setw(15) << setprecision (6) << Sample.meanE;
  // cout << setw(15) << setprecision (6) << Sample.meanE2;
  // cout << setw(15) << setprecision (6) << Sample.meanM;
  // cout << setw(15) << setprecision (6) << Sample.meanM2;
  // cout << setw(15) << setprecision (6) << Sample.absM << endl;
}

void printMatrix(sample Sample, int dim){
  cout << "==========================================================" << endl;
   for(int i = 0; i < dim+2; i++){
    for(int j = 0; j < dim+2; j++){
      cout << setw(12) << setprecision (10) << Sample.lattice[i][j];
    }
    cout << endl;
  }
  cout << "==========================================================" << endl;
}



/*
for(int t = 0; t < nTemps +1; t++){

  double dT = (T1 - T0)/(double) nTemps;
  double T = T0 + t*dT;

  sample S = sample(cycles, dim, T);
  S.initializeLattice(ordered);
  S.precalculateProbs();
  S.metropolis();
  S.normalize();
  printEnergies(S);
}
*/

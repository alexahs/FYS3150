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
// #include "initializeLattice.h"
// #include "analytical.h"
#include "sample.h"


using namespace std;

void output(sample Sample);

int main(){
  int dim = 2;
  int ordered = 1;
  double T0 = 2.0;
  double T1 = 3.0;
  int nTemps = 10;
  int cycles = 1<<20;

  for(int t = 0; t < nTemps +1; t++){
    double dT = (T1 - T0)/(double) nTemps;
    double T = T0 + t*dT;

    sample S = sample(cycles, dim, T);
    S.initializeLattice(ordered);
    S.precalculateProbs();
    S.metropolis();
    S.normalize();
    output(S);
  }

  return 0;
}

void output(sample Sample){
  cout << setw(2) << setprecision (6) << Sample.temperature;
  cout << setw(15) << setprecision (6) << Sample.meanE;
  cout << setw(15) << setprecision (6) << Sample.meanE2;
  cout << setw(15) << setprecision (6) << Sample.meanM;
  cout << setw(15) << setprecision (6) << Sample.meanM2;
  cout << setw(15) << setprecision (6) << Sample.absM << endl;
}

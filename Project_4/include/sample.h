#ifndef SAMPLE_H
#define SAMPLE_H
// #include <iostream>
// #include <cmath>
// #include <fstream>
// #include <iomanip>
// #include <cstdio>
// #include <chrono>
// #include <algorithm>
// #include <vector>
// #include <string>
// #include <time.h>
// #include <stdlib.h>
// #include <cmath>

// using namespace std;

class sample{


public:

  //Properties
  int dim;
  int nCycles;
  double temperature;
  double meanE;
  double meanE2;
  double meanM;
  double meanM2;
  double absM;
  double w[17];
  double **lattice;

  //Initializers
  // sample();
  sample(int cycs, int dim, double T);

  //Methods
  void initializeLattice(int ordered);
  void precalculateProbs();
  void metropolis();
  void normalize();
  double susceptibility();
  double heatCapacity();


};




#endif //SAMPLE_H

#ifndef SOLVER_H
#define SOLVER_H
#include "heavenlyBody.h"
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;



class solver{

public:
  //Properties
  int gridPoints;
  double simulationTime;
  vector <heavenlyBody> bodies;
  int solverType;

  //Initializers
  solver();
  solver(int n, double T, vector <heavenlyBody> planets, int m);

  //Methods
  void solve();

};



#endif // SOLVER_H

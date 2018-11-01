#include "analytical.h"
using namespace std;



double energyExpected(){
  double Z = 12.0 + 4.0*cosh(8.0);
  return -32.0/Z*sinh(8.0);
}

double energySquaredExpected(){
  double Z = 12 + 4*cosh(8.0);
  return 256.0/Z*cosh(8.0);
}

double magnExpected(){
  return 0.0;
}

double magnSquaredExpected(){
  double Z = 12 + 4*cosh(8.0);
  return 64.0/Z*(1 + exp(8.0));
}

double heatCapacity(){
  double Z = 12 + 4*cosh(8.0);
  256.0/(Z*Z)*(Z*cosh(8) - 2*sinh(8)*sinh(8));
}

double susceptibility(){
  double Z = 12 + 4*cosh(8.0);
  return 32.0/Z*(1 + exp(8));
}

#include <random>
#include <iostream>
#include <iomanip>
using namespace std;

double E_() {
    double Z = 4*cosh( 8 ) + 12;
    return -32/Z*sinh( 8 );
}
double E2_() {
    double Z = 4*cosh( 8 ) + 12;
    return 256/Z*cosh( 8 );
}
double M_() {
    return 0;
}
double M2_() {
    double Z = 4*cosh( 8 ) + 12;
    return 32/Z*(exp( 8 ) + 1 );
}
double chi() {
    return M2_()/* - M_()*M_()*/;
}
double Cv() {
    return E2_() - E_()*E_();
}
double absM() {
    double Z = 4*cosh( 8 ) + 12;
    return 8/Z*(exp( 8 ) + 2 );
}

void printAnalytical(){
  double d = 4.0;
  cout << setw(15) << setprecision(8) << "E";
  cout << setw(15) << setprecision(8) << "E2";
  cout << setw(15) << setprecision(8) << "M";
  cout << setw(15) << setprecision(8) << "M2";
  cout << setw(15) << setprecision(8) << "abs M";
  cout << setw(15) << setprecision(8) << "Cv";
  cout << setw(15) << setprecision(8) << "chi" << endl;
  cout << setw(15) << setprecision(8) << E_()/d;
  cout << setw(15) << setprecision(8) << E2_()/d;
  cout << setw(15) << setprecision(8) << M_();
  cout << setw(15) << setprecision(8) << M2_();
  cout << setw(15) << setprecision(8) << absM();
  cout << setw(15) << setprecision(8) << Cv();
  cout << setw(15) << setprecision(8) << chi() << endl;
}

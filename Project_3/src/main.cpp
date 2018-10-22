#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <vector>
#include "heavenlyBody.h"
#include <string>

using namespace std;
ofstream outfile;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

double pi = M_PI;

int main(int argc, char* argv[])
{

  char *outfilename;
  if(argc < 4){
    cout << "Bad usage: " << argv[0] <<
    " Read output filename, simulation time and grid points on same line" << endl;
    exit(1);
  }

  outfilename = argv[1];
  double T_ = atof(argv[2]);
  int n = atoi(argv[3]);
  double h = T_/((double) n);


  vector<double> accelerationOld(3);
  vector<double> accelerationNew(3);
  vector<heavenlyBody> allBodies;

  //mercury precession vectors
  // vector<double> peri_x;
  // vector<double> peri_y;
  // vector<double> peri_r;

  double FourPi2 = 4*pi*pi;
  double solarMass = 1.9891e30;

  // heavenlyBody Sun = heavenlyBody(1.0,0,0,0,0,0,0);
  // heavenlyBody Mercury = heavenlyBody(1.307e22/solarMass,-3.294363957786441E-01,-2.928799526088138E-01,5.618346324205380E-03,365.25*1.320405892727915E-02, 365.25*-1.952252048338632E-02, 365.25*-2.807294373094382E-03);
  // heavenlyBody Venus = heavenlyBody(4.867E24/solarMass, 7.243545955158947E-01, -3.278712379892032E-02, -4.242728890559555E-02, 365.25*1.017391327967621E-03, 365.25* 2.010584861519629E-02, 365.25* 2.168289888508737E-04);
  // heavenlyBody Earth = heavenlyBody(5.972E24/solarMass,9.837576984919719E-01,1.889233809711713E-01,-8.631011464030984E-05,-3.406523137555859E-03*365.25,1.686035619678342E-02*365.25,-1.194254105980157E-06*365.25);
  // heavenlyBody Mars = heavenlyBody(6.39E23/solarMass, 1.349004548617385E+00, -2.975589233646888E-01, -3.956440841859040E-02, 365.25*3.610034965148588E-03, 365.25* 1.484808760059448E-02, 365.25* 2.224945616221949E-04);
  // heavenlyBody Jupiter = heavenlyBody(1.898E27/solarMass,-2.724865762194714E+00,-4.624789318060123E+00,8.013249787610907E-02,6.411862928919486E-03*365.25,-3.471082490961821E-03*365.25,-1.290147901227175E-04*365.25);
  // heavenlyBody Saturn = heavenlyBody(5.683E26/solarMass, 1.497082568921199E+00, -9.943470921581483E+00,  1.132983557425057E-01, 365.25*5.209583578051823E-03, 365.25* 8.120803848912152E-04, 365.25*-2.211308505468577E-04);
  // heavenlyBody Uranus = heavenlyBody(8.681E25/solarMass, 1.719595695177778E+01,  9.965486713193039E+00, -1.857636424997038E-01, 365.25*-2.000761535443054E-03, 365.25* 3.219594226509228E-03, 365.25* 3.795711294500967E-05);
  // heavenlyBody Neptune = heavenlyBody(1.024E26/solarMass, 2.891239407445751E+01, -7.753050308782163E+00, -5.066556247342422E-01, 365.25*7.926104454699854E-04, 365.25* 3.050689379330089E-03, 365.25*-8.139915196891708E-05);
  // heavenlyBody Pluto = heavenlyBody(1.30900E22/solarMass, 1.161374129179143E+01, -3.157937303069106E+01,  1.979427629835602E-02, 365.25*3.006977217402132E-03, 365.25* 4.205759240708480E-04, 365.25*-9.057561756443009E-04);
  // allBodies.push_back(Sun);
  // allBodies.push_back(Mercury);
  // allBodies.push_back(Venus);
  // allBodies.push_back(Earth);
  // allBodies.push_back(Mars);
  // allBodies.push_back(Jupiter);
  // allBodies.push_back(Saturn);
  // allBodies.push_back(Uranus);
  // allBodies.push_back(Neptune);
  // allBodies.push_back(Pluto);


  //mercury precession temporary variables
  // double temp0 = 0;
  // double temp1;
  // double temp2;
  // double tempx;
  // double tempy;

  outfile.open(outfilename);
  for(int i = 0; i < n; i++)
  {
    for(heavenlyBody &current: allBodies)
    {

        //Mercury precession calculations
        // if(current.mass != 1.0){
        //   temp1 = current.distance(sun);
        // }

        //Set acceleration to 0
        accelerationOld = {0.0, 0.0, 0.0};
        accelerationNew = {0.0, 0.0, 0.0};

        //Calculate acceleration
        for(heavenlyBody &other: allBodies)
        {
          accelerationOld[0] += current.newtonian_acceleration(other, 0);
          accelerationOld[1] += current.newtonian_acceleration(other, 1);
          accelerationOld[2] += current.newtonian_acceleration(other, 2);
        }

        //Calculate new position
        current.position[0] = current.position[0] + h*current.velocity[0] + 0.5*h*h*accelerationOld[0];
        current.position[1] = current.position[1] + h*current.velocity[1] + 0.5*h*h*accelerationOld[1];
        current.position[2] = current.position[2] + h*current.velocity[2] + 0.5*h*h*accelerationOld[2];


        //Mercury precession calculations
        // if(current.mass != 1.0){
        //   temp2 = current.distance(sun);
        //   tempx = current.position[0];
        //   tempy = current.position[1];
        // }
        //
        // if(current.mass != 1.0 && temp0 > temp1 && temp2 > temp1){
        //   peri_x.push_back(tempx);
        //   peri_y.push_back(tempy);
        // }
        //
        // if(current.mass != 1.0){
        //   temp0 = temp1;
        // }


        //Sun fixed in the origin
        // allBodies[0].position[0] = 0.0;
        // allBodies[0].position[1] = 0.0;
        // allBodies[0].position[2] = 0.0;

        //Calculate acceleration in next step
        for(heavenlyBody &other: allBodies)
        {
          accelerationNew[0] += current.newtonian_acceleration(other, 0);
          accelerationNew[1] += current.newtonian_acceleration(other, 1);
          accelerationNew[2] += current.newtonian_acceleration(other, 2);
        }

        //Calculate new velocity
        current.velocity[0] = current.velocity[0] + 0.5*h*(accelerationOld[0] + accelerationNew[0]);
        current.velocity[1] = current.velocity[1] + 0.5*h*(accelerationOld[1] + accelerationNew[1]);
        current.velocity[2] = current.velocity[2] + 0.5*h*(accelerationOld[2] + accelerationNew[2]);

        //Write to file
        outfile << setw(30) << setprecision(12) << current.position[0];
        outfile << setw(30) << setprecision(12) << current.position[1];
        // outfile << setw(30) << setprecision(12) << current.position[2];
  }
    outfile << endl;
  }
  outfile.close();


  //Mercury precession file writing
  // int array_length = 0;
  // for(double &current: peri_x){
  //   array_length += 1;
  // }
  // cout << array_length << endl;
  // outfile.open("precession_newtonian.dat");
  // for(int i = 0; i < array_length; i++){
  //   outfile << setw(20) << setprecision(12) << peri_x[i];
  //   outfile << setw(20) << setprecision(12) << peri_y[i] << endl;
  // }
  // outfile.close();



}

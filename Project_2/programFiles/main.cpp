#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <time.h>
#include <chrono>
#include <armadillo>
#include "methods.h"
#include "unitTests.h"
using namespace std;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
ofstream ofile;

void solver(int n, int electrons, double rhoMax, double omega);



int main(int argc, char *argv[]){
    //Commandline arguments needed to run this program are as follows:
    //1: number of grid points, 2: specifies which problem to solve, 0, 1 or 2
    //3: specifies rhoMax (for quantum dots problems)
    //4: specifies strength of harmonic oscillator in the two electron problem
    if(argc < 4){
      cout << "Four commandline arguments are needed." << endl;
      exit(1);
    }

    int n = atoi(argv[1]);
    int electrons = atoi(argv[2]);
    double rhoMax = atof(argv[3]);
    double omega = atof(argv[4]);

    solver(n, electrons, rhoMax, omega);

    // testMethods(n);

    return 0;

}



void solver(int n, int electrons, double rhoMax, double omega){
    //This function sets up values for an upper triangular matrix A corresponding
    //to the two electron problem and uses the methods offdiag and jacobi_rotate
    //for finding its eigenvalues. The variable electrons determines what kind of problem
    //we are looking at, with 0 being the buckling beam problem.


    //Initialize matrices and constants
    double **A = new double*[n];
    double **R = new double*[n];
    double *rho = new double[n];
    double rho0 = 0.0;
    double *d = new double[n];
    double h = (rhoMax - rho0)/((double) n - 1);
    double a;
    double *eigenvalues = new double[n];

    for(int i = 0; i < n; i++){
      A[i] = new double[n];
      R[i] = new double[n];
    }


    //Calculate values for the diagonals of matrix A
    if(electrons == 2){
      a = -1.0/(h*h);
      for(int i = 0; i < n; i++){
        rho[i] = rho0 + (i+1)*h;
        d[i] = 2.0/(h*h) + omega*omega*rho[i]*rho[i] + 1.0/rho[i];
      }
    }
    else if(electrons == 1){
      a = -1.0/(h*h);
      for(int i = 0; i < n; i++){
        rho[i] = rho0 + (i+1)*h;
        d[i] = 2.0/(h*h) + rho[i]*rho[i];
      }
    }
    else{
      a = -1.0/(h*h);
      for(int i = 0; i < n; i++){
        rho[i] = rho0 + (i+1)*h;
        d[i] = 2.0/(h*h);
      }
    }




    //Set up initial matrix values
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A[i][j] = 0.0;
            R[i][j] = 0.0;
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j){
                A[i][j] = d[i];
                R[i][j] = 1.0;
            }
            else if (i == j-1 || i == j+1){
                A[i][j] = a;
              }
        }
    }


    //Tolerance values and iteration loop
    double tolerance = 1.0E-10;
    int iterations = 0;
    int maxiterations = n*n*n;
    double maxval = 1000.0;

    steady_clock::time_point programStart;
    programStart = steady_clock::now();
    while(maxval > tolerance && iterations < maxiterations){
        int p, q;
        offdiag(A, n, p, q, &maxval);
        jacobi_rotate(A, R, p, q, n);
        iterations++;
    }
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    cout << "Program time: " << programTime.count() << " seconds" << endl;
    cout << "no. of jacobi iterations: " << iterations << endl;

    //Store eigenvalues in a vector
    for(int i = 0; i < n; i++){
      eigenvalues[i] = A[i][i];
    }

    //find index of lowest eigenvalue
    double min_element = 1000.0;
    int k;
    for(int i = 0; i < n; i++){
      if(eigenvalues[i] < min_element){
        min_element = eigenvalues[i];
        k = i;
      }
    }

    sort(eigenvalues, eigenvalues + n);

    cout << "eigenvalues:" << endl;
    for(int i = 0; i < 4; i++){
      cout << setw(12) << setprecision (10) << eigenvalues[i] << endl;
    }


    //Deallocate memory
    delete[] rho;
    delete[] d;
    delete[] eigenvalues;
    for(int i = 0; i < n; i++){
      delete[] A[i];
      delete[] R[i];
    }
    delete[] A;
    delete[] R;




}

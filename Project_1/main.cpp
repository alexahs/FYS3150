#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include "time.h"
#include <lib.h>
#include <chrono>

using namespace std;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;



ofstream ofile;

double generalAlgo(int n);
double specializedAlgo(int n);
double maxError(int n);
double LUdecomp(int n);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);

int main(){



    int n = 1000;
    int N = 10;

//    int N = 4;
    specializedAlgo(n);
//    generalAlgo(n);
//    maxError(n);
//    LUdecomp(n);




    return 0;
}

double specializedAlgo(int n){
    double h = 1/(double(n));

    // Initialize vectors
    double *bVec = new double[n+1];
    double *bTilde = new double[n+1];
    double *fTilde = new double[n+1];
    double *vVec = new double[n+1];
    double *fVec = new double[n+1];
    double *x = new double[n+1];
    double *exact = new double[n+1];

    // Set boundry conditions
    for(int i = 0; i < n+1; i++){
        double j = i;
        fVec[i] = h*h*100*exp(-10*(j*h));
        x[i] = i*h;

    }

    for(int i = 0; i < n; i++){
        bVec[i] = 2.;
    }


    bTilde[1] = 2.0;
    fTilde[1] = fVec[1];
    vVec[0] = 0;
    vVec[n] = 0;
    vVec[n] = 0;

    // Analytical solution
    for(int i = 0; i < n+1; i++){
        exact[i] = 1 - (1 - exp(-10))*x[i] - exp(-10*x[i]);
    }



    steady_clock::time_point programStart;
    programStart = steady_clock::now();

    // Forward substitution
    for(int i = 2; i < n+1; i++){
        bTilde[i] = bVec[i] - 1.0/bTilde[i-1];
        fTilde[i] = fVec[i] + fTilde[i-1]/bTilde[i-1];

    }
    //Backward substitution
    for(int i = n-1; i > 0; i--){
        vVec[i] = (fTilde[i] + vVec[i+1])/bTilde[i];
    }

    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);




    // Write results to file
    ofile.open("specialized_results_n1000.txt");
    for(int i = 0; i < n + 1; i++){
        ofile << setw(8) << x[i];
        ofile << setw(20) << setprecision(8) << exact[i];
        ofile << setw(20) << setprecision(8) << vVec[i];
//        ofile << setw(20) << setprecision(8) << Q << endl;
    }
    ofile.close();

    delete [] bVec;
    delete [] fTilde;
    delete [] vVec;
    delete [] fVec;
    delete [] exact;

    return programTime.count();

}

double maxError(int n){

    double *errors = new double[n+1];
    double *stepSizes = new double[n+1];

    double LUd = LUdecomp(n);

    double emax = 0.0;

//    for(int i = 0; i < n+1; i++){
//        if(LUd[i] > emax){
//            emax = LUd[i];
//        }
//    }


//    emax = log10(emax);
//    cout << "n=" << n << "   emax=" << emax << endl;


    delete[] errors;
    delete[] stepSizes;
}

double generalAlgo(int n){
    double h = 1/(double(n));

    // Initialize vectors
    double *aVec = new double[n+1];
    double *bVec = new double[n+1];
    double *bTilde = new double[n+1];
    double *fTilde = new double[n+1];
    double *cVec = new double[n+1];
    double *vVec = new double[n+1];
    double *fVec = new double[n+1];
    double *x = new double[n+1];
    double *exact = new double[n+1];
    double *relErrors = new double[n+1];
    double K;

    for(int i = 0; i < n+1; i++){
        double j = i;
        fVec[i] = h*h*100*exp(-10*(j*h));
        x[i] = i*h;

    }

    for(int i = 0; i < n+1; i++){
        aVec[i] = -1.;
        cVec[i] = -1.;
        bVec[i] = 2.;
    }

    // Set boundry conditions
    bTilde[1] = bVec[1];
    fTilde[1] = fVec[1];
    vVec[0] = 0;
    vVec[n] = 0;
    vVec[n-1] = fTilde[n-1]/bTilde[n-1];

    // Analytical solution
    for(int i = 0; i < n; i++){
        exact[i] = 1 - (1 - exp(-10))*x[i] - exp(-10*x[i]);
    }

    steady_clock::time_point programStart;
    programStart = steady_clock::now();
    // Forward substitution
    for(int i = 2; i < n+1; i++){
        K = aVec[i-1]/bTilde[i-1];
        bTilde[i] = bVec[i] - cVec[i-1]*K;
        fTilde[i] = fVec[i] - fTilde[i-1]*K;
    }

    //Backward substitution
    for(int i = n-1; i > 0; i--){
        vVec[i] = (fTilde[i] - cVec[i]*vVec[i+1])/bTilde[i];
    }

    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);





    //Calculate error
    for(int i = 0; i < n+1; i++){
        relErrors[i] = log10(fabs((exact[i] - vVec[i])/exact[i]));
    }

    //Find max error
    double emax = 0.0;
    for(int i = 1; i < n; i++){
        if(relErrors[i] < emax){
            emax = relErrors[i];
        }
    }

    // Write results to file
//    ofile.open("generalized_results_n1000.txt");
//    for(int i = 0; i < n+1 ; i++){
//        ofile << setw(8) << x[i];
//        ofile << setw(25) << setprecision(15) << exact[i];
//        ofile << setw(25) << setprecision(15) << vVec[i];
//        ofile << setw(25) << setprecision(15) << relErrors[i] << endl;

//    }
//    ofile.close();



    delete [] aVec;
    delete [] bVec;
    delete [] bTilde;
    delete [] fTilde;
    delete [] cVec;
    delete [] vVec;
    delete [] fVec;
    delete [] x;
    delete [] exact;
    delete [] relErrors;

    return programTime.count();

}

void ludcmp(double **a, int n, int *indx, double *d){
   int      i, imax, j, k;
   double   big, dum, sum, temp, *vv;

  vv = new(nothrow) double [n];
  if(!vv) {
    printf("\n\nError in function ludcm():");
    printf("\nNot enough memory for vv[%d]\n",n);
    exit(1);
  }

   *d = 1.0;                              // no row interchange yet
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = ZERO;
      for(j = 0; j < n; j++) {
         if((temp = fabs(a[i][j])) > big) big = temp;
      }
      if(big == ZERO) {
         printf("\n\nSingular matrix in routine ludcmp()\n");
         exit(1);
      }
      vv[i] = 1.0/big;                 // save scaling */
   } // end i-loop */

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i< j; i++) {   // not i = j
         sum = a[i][j];
     for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
     a[i][j] = sum;
      }
      big = ZERO;   // initialization for search for largest pivot element
      for(i = j; i< n; i++) {
         sum = a[i][j];
     for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
     a[i][j] = sum;
     if((dum = vv[i]*fabs(sum)) >= big) {
        big = dum;
        imax = i;
     }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0;k< n; k++) {       // yes
        dum        = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k]    = dum;
     }
     (*d)    *= -1;            // and change the parit of d
     vv[imax] = vv[j];         // also interchange scaling factor
      }
      indx[j] = imax;
      if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;

        /*
        ** if the pivot element is zero the matrix is singular
        ** (at least to the precision of the algorithm). For
        ** some application of singular matrices, it is desirable
        ** to substitute ZERO for zero,
        */

      if(j < (n - 1)) {                   // divide by pivot element
         dum = 1.0/a[j][j];
     for(i=j+1;i < n; i++) a[i][j] *= dum;
      }
   } // end j-loop over columns

   delete [] vv;   // release local memory
} // End: function ludcmp()

void lubksb(double **a, int n, int *indx, double *b){
   int        i, ii = -1, ip, j;
   double     sum;

   for(i = 0; i< n; i++) {
      ip    = indx[i];
      sum   = b[ip];
      b[ip] = b[i];
      if(ii > -1)   for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
      else if(sum) ii = i;
      b[i] = sum;
   }
   for(i = n - 1; i >= 0; i--) {
      sum = b[i];
      for(j = i+1; j < n; j++) sum -= a[i][j] * b[j];
      b[i] = sum/a[i][i];
   }
} // End: function lubksb()

double LUdecomp(int n){


//    cout << "Initializing vectors.." << endl;
    //Initialize matrix and vectors
    double **A = new double*[n];
    int *indices = new int[n];
    double *f = new double[n];
    double *x = new double[n];
    double h = 1/(double(n));
    double *d = new double[n];
    double *anal = new double[n];
    double *relErrors = new double[n];

    for(int i = 0; i < n; i++){
        A[i] = new double[n];
    }

    //Values for steps x and known function f
    for(int i = 0; i < n+1; i++){
        double j = i;
        f[i] = h*h*100*exp(-10*(j*h));
        x[i] = i*h;
    }

    // Analytical solution
    for(int i = 0; i < n; i++){
        anal[i] = 1 - (1 - exp(-10))*x[i] - exp(-10*x[i]);
    }

//    cout << "Setting up matrix values.." << endl;

    //Set up values for matrix
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A[i][j] = 0.0;
        }
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j){A[i][j] = 2.0;}
            else if (i == j-1 || i == j+1){A[i][j] = -1.0;}
        }
    }

//    cout << "Computing LU decomposition.." << endl;

    steady_clock::time_point programStart;
    programStart = steady_clock::now();

    ludcmp(A, n, indices, d);

//    cout << "Computing back substitution.." << endl;


    lubksb(A, n, indices, f);


    f[0] = 0.0;
    f[n] = 0.0;

    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);
    //Calculate error
    for(int i = 0; i < n+1; i++){

        relErrors[i] = fabs((anal[i] - f[i])/anal[i]);
//        cout << relErrors[i] << endl;
    }

    cout << "time...." << programTime.count() << endl;
    return programTime.count();






    // Write results to file
    ofile.open("LUdecomp_results_n1000.txt");
    for(int i = 0; i < n+1 ; i++){
        ofile << setw(8) << setprecision(15) << x[i];
        ofile << setw(25) << setprecision(15) << anal[i];
        ofile << setw(25) << setprecision(15) << f[i];
        ofile << setw(25) << setprecision(15) << relErrors[i] << endl;


    }
    ofile.close();



}

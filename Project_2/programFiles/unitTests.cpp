#include "methods.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <chrono>
#include <algorithm>

using namespace std;

void testJacobi_rotate(int n);
void testOffdiag(int n);


void testMethods(int n){

  testOffdiag(n);
  testJacobi_rotate(n);

}


void testJacobi_rotate(int n){
    //this function test the jaboci rotate method of diagonalizing a matrix and finding its eigenvalues.
    //In this test case the analytical solution of the eigenvalues are known, and are compared to the
    //computed ones to determine if the function is working properly.


    //Initialize matrices and vectors
    double pi = acos(-1.0);
    double **A = new double*[n];
    double **R = new double*[n];
    double *jacobiEigenvals = new double[n];
    double *analyticalEigenvals = new double[n];



    for(int i = 0; i < n; i++){
      A[i] = new double[n];
      R[i] = new double[n];
    }

    //Analytical eigenvalues
    for(int i = 1; i < n+1; i++){
      analyticalEigenvals[i-1] = 2 + -2*cos(i*pi/(n+1.0));
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
                A[i][j] = 2.0;
                R[i][j] = 1.0;
            }
            else if (i == j-1 || i == j+1){
                A[i][j] = -1.0;
              }
        }
    }


    //rotates the matrix untill all non-diagonal elements are 0.
    double tolerance = 1.0E-10;
    int iterations = 0;
    int maxiterations = n*n*n;
    double maxval = 1000.0;
    while(maxval > tolerance && iterations < maxiterations){
        int p, q;
        offdiag(A, n, p, q, &maxval);
        jacobi_rotate(A, R, p, q, n);
        iterations++;
    }


    //stores the eigenvalue in a single vector.
    for(int i = 0; i < n; i++){
      jacobiEigenvals[i] = A[i][i];
    }

    sort(jacobiEigenvals, jacobiEigenvals + n);
    sort(analyticalEigenvals, analyticalEigenvals + n);


    double eps = 1e-5;
    double sumError = 0.0;


    //compares analytical eigenvalues to the computed eigenvalues.
    for(int i = 0; i < n; i++){
      sumError += fabs(jacobiEigenvals[i] - analyticalEigenvals[i]);
    }


    if(sumError > eps){
      cout << "Error in jacobi_rotate method." << endl;
    }
    else{
      cout << "testJacobi_rotate complete, working as intended." << endl;
    }


}

void testOffdiag(int n){
    //This function tests the method for finding the matrix with the biggest value.
    //It creates a matrix with known indices of the biggest matrix element and compares
    //these to the indices returned by the offdiag function.



    //predetermined indices for the biggest matrix element
    int index1 = 5;
    int index2 = 6;
    int P, Q;


    //initialize matrix
    double **A = new double*[n];
    for(int i = 0; i < n; i++){
      A[i] = new double[n];
    }
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        A[i][j] = i+j;
      }
    }

    //value of larges matrix element
    A[index1][index2] = 1000.0;

    double tolerance = 1.0E-10;
    int iterations = 0;
    int maxiterations = n*n*n;
    double maxval = 1000.0;
    while(maxval > tolerance && iterations < maxiterations){
        int p, q;
        P = p;
        Q = q;
        offdiag(A, n, p, q, &maxval);
        iterations++;
      }


    if(index1 == P && index2 == Q){
        cout << "testOffdiag complete, working as intended." << endl;
    }
    else{
        cout << "Error in testOffdiag." << endl;

    }

}

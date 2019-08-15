#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <chrono>
#include <armadillo>
#include <algorithm>


using namespace arma;
using namespace std;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;



void offdiag(double **A, int n, int &p, int &q, double *maxval);
void jacobi_rotate(double **A, double **R, int p, int q, int n);



void offdiag(double **A, int n, int &p, int &q, double *maxval){
    //This method searches though an upper triangular matrix A of dimensions
    //nxn and returns the indices of the largest element in A. It returns both
    //the indices of this element, aswell as the value of the element itself.
    double max = 0.0;
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; ++j){
            double aij = fabs(A[i][j]);
            if(aij > max){
                max = aij;
                p = i;
                q = j;
            }
        }
    }
    *maxval = max;

}

void jacobi_rotate(double **A, double **R, int k, int l, int n){

    //This is a method for transforming an upper triangular matrix A
    //to a diagonal matrix. This is done by chosing the largest matrix element A[k][l]
    //and rotating the matrix an angle theta, such that A[k][l] becomes zero.
    //The indices k and l are provided by the offdiag method.
    //R is the matrix for storing the eigenvectors.

    double s, c, t, tau;
    if(A[k][l] != 0.0){
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);
        if(tau >= 0){
          t = 1.0/(tau + sqrt(1 + tau*tau));
        }
        else{
            t = -1.0/(-tau + sqrt(1 + tau*tau));
        }
        c = 1.0/(sqrt(1 + t*t));
        s = t*c;
    }
    else{
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    A[k][k] = c*c*a_kk - 2*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = c*c*a_ll + 2*c*s*A[k][l] + s*s*a_kk;
    A[k][l] = 0.0;
    A[l][k] = 0.0;
    for(int i = 0; i < n; i++){
        if(i != k && i != l){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[i][l] = c*a_il + s*a_ik;
            A[k][i] = A[i][k];
            A[l][i] = A[i][l];
        }
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
}

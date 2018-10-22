#ifndef METHODS_H
#define METHODS_H


void offdiag(double **A, int n, int &p, int &q, double *maxval);
void jacobi_rotate(double **A, double **R, int p, int q, int n);


#endif

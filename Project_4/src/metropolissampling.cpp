#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "initializelattice.h"
using namespace std;

void MetropolisSampling ( int dim, int MCcycles, int loopStart, int loopStop, double T, double *ExpectVal, int ordered) {
    //  Random number
    random_device rd;
    mt19937_64 gen( rd() );
    //  Uniform "REAL" distribution for x \in [0,1]
    uniform_real_distribution<double>RandomNumberGenerator( 0, 1 );
    uniform_int_distribution<int>SpinDistribution( 0, dim-1 );


    double **SpinMatrix = new double*[dim];
    for ( int i = 0; i < dim; i++ ) SpinMatrix[i] = new double[dim];
    double E = 0;
    double M = 0;
    for ( int i = 0; i < 5; i++ ) ExpectVal[i] = 0;

    InitializeLattice ( dim, SpinMatrix, E, M , ordered);
    double *EnergyDifference = new double[17];

    for ( int i = 0; i < 17; i += 4 ) EnergyDifference[i] = exp( -( i-8 )/T );

    //  Start Monte Carlo cycle
    for ( int cycle = loopStart; cycle <= loopStop; cycle++ ) {
        for ( int x = 0; x < dim; x++ ) {
            for ( int y = 0; y < dim; y++ ) {
                int ix = SpinDistribution( gen );
                int iy = SpinDistribution( gen );
                int DeltaE = 2*SpinMatrix[ix][iy]*
                        ( SpinMatrix[ix][PeriodicBoundary( iy, dim, -1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, -1 )][iy] +
                        SpinMatrix[ix][PeriodicBoundary( iy, dim, 1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, 1 )][iy]);
                if ( RandomNumberGenerator( gen ) <= EnergyDifference[DeltaE + 8] ) {
                    SpinMatrix[ix][iy] *= -1.0; // Flip one spin => new configuration
                    M += 2*SpinMatrix[ix][iy];
                    E += DeltaE;
                }
            }
        }
        // cout << ExpectVal[0] << endl;
        ExpectVal[0] += E;
        ExpectVal[1] += E*E;
        ExpectVal[2] += M;
        ExpectVal[3] += M*M;
        ExpectVal[4] += fabs( M );

    }
    double norm = 1/(( double ) MCcycles );
    double meanE  = ExpectVal[0]*norm;
    double meanE2 = ExpectVal[1]*norm;
    double meanM  = ExpectVal[2]*norm;
    double meanM2 = ExpectVal[3]*norm;
    double absM   = ExpectVal[4]*norm;
    double varE   = ( meanE2 - meanE*meanE )/( dim*dim );
    double varM   = ( meanM2 - meanM*meanM )/( dim*dim );

//    cout << endl;
//    cout << "L                = " << dim << endl;
//    cout << "T                = " << T << endl;
//    cout << "Number of cycles = " << MCcycles << endl << endl;

//    cout << "Metropolis gives " << endl;
//    cout << "<E>              = " << meanE << endl;
//    cout << "|M|              = " << absM << endl;
//    cout << "Heat capacity    = " << varE << endl;
//    cout << "Susceptibility   = " << varM << endl;

//    ofstream outfile;
//    outfile.open("info.txt");
//    outfile << setw(15) << setprecision(8) << T;
//    outfile << setw(15) << setprecision(8) << meanE;
//    outfile << setw(15) << setprecision(8) << meanE2;
//    outfile << setw(15) << setprecision(8) << meanM;
//    outfile << setw(15) << setprecision(8) << meanM2;
//    outfile << setw(15) << setprecision(8) << absM;
//    outfile << setw(15) << setprecision(8) << varM/T;
//    outfile << setw(15) << setprecision(8) << varE/T/T;
//    outfile.close();

    //  Memory deallocation
    for ( int i = 0; i < dim; i++ ) {
        delete [] SpinMatrix[i];
    }
    delete [] SpinMatrix;
    delete [] EnergyDifference;
}

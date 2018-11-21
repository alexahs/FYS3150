
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "initializelattice.h"
using namespace std;
ofstream ofile;

void MetropolisSampling ( int dim, int MCcycles, int loopStart, int loopStop, double T, double *ExpectVal,
                          int ordered, double &Accept ) {
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

    InitializeLattice ( dim, SpinMatrix, E, M, ordered );
    double *EnergyDifference = new double[17];

    //Initialize delta E vector
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
                        SpinMatrix[PeriodicBoundary( ix, dim, 1 )][iy] );
                if ( RandomNumberGenerator( gen ) <= EnergyDifference[DeltaE + 8] )
                    //accept new configuration
                    SpinMatrix[ix][iy] *= -1.0;
                    M += 2*SpinMatrix[ix][iy];
                    E += DeltaE;
                }
            }
        }//end lattice sweep

        //Update expectation values
        ExpectVal[0] += E;
        ExpectVal[1] += E*E;
        ExpectVal[2] += M;
        ExpectVal[3] += M*M;
        ExpectVal[4] += fabs( M );

    } //end MC loop


    //  Memory deallocation
    for ( int i = 0; i < dim; i++ ) {
        delete [] SpinMatrix[i];
    }
    delete [] SpinMatrix;
    delete [] EnergyDifference;
}

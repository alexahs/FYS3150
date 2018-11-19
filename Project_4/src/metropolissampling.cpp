
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "initializelattice.h"
using namespace std;
ofstream ofile;

void MetropolisSampling ( int dim, int MCcycles, int loopStart, int loopStop, double T, double *ExpectVal,
                          int ordered, double &Accept ) {
    //  Random number
    random_device rd;
    mt19937_64 gen( rd() );
    //  Uniform "REAL" distribution for x \in [0,1]
    uniform_real_distribution<double>RandomNumberGenerator( 0, 1 );
    //  Initialize lattice spin values
    uniform_int_distribution<int>SpinDistribution( 0, dim-1 );

    double **SpinMatrix = new double*[dim];
    for ( int i = 0; i < dim; i++ ) SpinMatrix[i] = new double[dim];
    double E = 0;
    double M = 0;
    for ( int i = 0; i < 5; i++ ) ExpectVal[i] = 0;

    InitializeLattice ( dim, SpinMatrix, E, M, ordered );
    double *EnergyDifference = new double[17];

    //  4c --
    vector<vector<int>> EMV;
    EMV.resize(2);
    int numOfAccepts = 0;
    // EMV[0].push_back( (int) E );
    // EMV[1].push_back( (int) M );
    //  --
    for ( int i = 0; i < 17; i += 4 ) EnergyDifference[i] = exp( -( i-8 )/T );
    //  Start Monte Carlo cycle
    for ( int cycle = loopStart; cycle <= loopStop; cycle++ ) {
        for ( int x = 0; x < dim; x++ ) {
            for ( int y = 0; y < dim; y++ ) {
                int ix = SpinDistribution( gen );
                int iy = SpinDistribution( gen );
                //  Edge elements are 'ignored' using the next couple of lines
                int DeltaE = 2*SpinMatrix[ix][iy]*
                        ( SpinMatrix[ix][PeriodicBoundary( iy, dim, -1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, -1 )][iy] +
                        SpinMatrix[ix][PeriodicBoundary( iy, dim, 1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, 1 )][iy] );
                if ( RandomNumberGenerator( gen ) <= EnergyDifference[DeltaE + 8] ) {
                    SpinMatrix[ix][iy] *= -1.0; // Flip one spin => new configuration
                    M += 2*SpinMatrix[ix][iy];
                    E += DeltaE;
                    numOfAccepts++;
                }
            }
        }
        // EMV[0].push_back( E ); //  push_back in c++ is like append i python
        // EMV[1].push_back( fabs( M ) );

        ExpectVal[0] += E;
        ExpectVal[1] += E*E;
        ExpectVal[2] += M;
        ExpectVal[3] += M*M;
        ExpectVal[4] += fabs( M );

    }
//    ofstream ofile;
    //  4c) write to binary file
    // ofile.open("energiesEachMCsample" + to_string(T) + "Cycles" + to_string(MCcycles) + ".bin", ofstream::binary);
    //
    // ofile.write(reinterpret_cast<const char*> (EMV[0].data()), EMV[0].size()*sizeof(int));
    // ofile.close();

    Accept = numOfAccepts/double(loopStop);

//    double norm = 1/(( double ) MCcycles );
//    double meanE  = ExpectVal[0]*norm;
//    double meanE2 = ExpectVal[1]*norm;
//    double meanM  = ExpectVal[2]*norm;
//    double meanM2 = ExpectVal[3]*norm;
//    double absM   = ExpectVal[4]*norm;
//    double varE   = ( meanE2 - meanE*meanE )/( dim*dim );
//    double varM   = ( meanM2 - meanM*meanM )/( dim*dim );

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

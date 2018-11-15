
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "analytical.h"
#include "metropolissampling.h"
#include <mpi.h>
#include <chrono>

using namespace std;
ofstream outfile;

void output( int dim, double T, double *ExpectVal, int MCcycles, double timing );

int main( int argc, char *argv[] )
{
    double *ExpectVal = new double[5];
    double *TotalExpectVal = new double[5];
    for( int i = 0; i < 5; i++ ) TotalExpectVal[i] = 0;

    //-------------------------------------------------------------------------
    //    Project 4c)
    string filename;
    int ordered = 1; //  Choose 1 for ordered matrix, choose 0 for random matrix
    int dim = 100;   //  Dimension of the matrix L
    int MCcycles = 1e5;
    double InitialTemp = 2.0;
    double FinalTemp = 2.5;
    double TimeStep = 0.05;
    double timing;
    chrono::high_resolution_clock::time_point t1;
    chrono::high_resolution_clock::time_point t2;


    //  Initialize parallellization
    int nProcs;
    int my_rank;

    MPI_Init (&argc, &argv);
    MPI_Comm_size ( MPI_COMM_WORLD, &nProcs );
    MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );


    //  Broadcast variables to all nodes on my CrapBook Air
    // MPI_Bcast (&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast (&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast (&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast (&TimeStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int cycleInterval = MCcycles/nProcs;
    int loopStart = my_rank*cycleInterval;
    int loopStop = (my_rank+1)*cycleInterval;


    if (my_rank == 0){
        t1 = chrono::high_resolution_clock::now();
        outfile.open("CriticalTemps_Dim" + to_string(dim) + "Cycles" + to_string(MCcycles) + ".dat", std::ios_base::app);
        outfile << setw(15) << setprecision(8) << "T";
        outfile << setw(15) << setprecision(8) << "E";
        outfile << setw(15) << setprecision(8) << "M";
        outfile << setw(15) << setprecision(8) << "C_V";
        outfile << setw(15) << setprecision(8) << "chi" << endl;
      }

    for ( double T = InitialTemp; T <= FinalTemp; T += TimeStep) {
        // cout << "HELLO MONEY \n";
        double acceptRatio = 0;

        if (my_rank==0) t1 = chrono::high_resolution_clock::now();
        MetropolisSampling( dim, MCcycles, loopStart, loopStop, T, ExpectVal, ordered, acceptRatio );
        for ( int i = 0; i < 5; i++ ) {
            MPI_Reduce(&ExpectVal[i], &TotalExpectVal[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        if (my_rank==0) {
            // timing = time_span.count();
            output( dim, T, TotalExpectVal, MCcycles, timing );
            cout << "T = " << T << " done...\n";
        }
//        ofstream outfile;
//        outfile.open("/Users/hennoz/FYS3150Project4/acceptsRatioVsT.txt", std::ios_base::app);
//        outfile << setw(15) << setprecision(8) << T;
//        outfile << setw(15) << setprecision(8) << acceptRatio << endl;
//        outfile.close();
    }
    outfile.close();
    if(my_rank == 0){
        t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = std::chrono::duration_cast<chrono::duration<double>>(t2 - t1);
        cout << "Program time: " << time_span.count() << endl;
    }

    delete [] ExpectVal;
    delete [] TotalExpectVal;
    MPI_Finalize ();
    return 0;
}


void output( int dim, double T, double *ExpectVal, int MCcycles, double timing ) {

  for( int i = 0; i < 5; i++ ) ExpectVal[i] /= (MCcycles);

  double E_variance = (ExpectVal[1] - ExpectVal[0]*ExpectVal[0])/dim/dim;
  double M_variance = (ExpectVal[3] - ExpectVal[4]*ExpectVal[4])/dim/dim;

  outfile << setw(15) << setprecision(8) << T;                   //Temperature
  outfile << setw(15) << setprecision(8) << ExpectVal[0]/dim/dim;        //Mean energy
  outfile << setw(15) << setprecision(8) << ExpectVal[4]/dim/dim;        //Mean abs magnetization
  outfile << setw(15) << setprecision(8) << E_variance/T/T;       //Heat capaticy
  outfile << setw(15) << setprecision(8) << M_variance/T << endl;   //Susceptibility
}

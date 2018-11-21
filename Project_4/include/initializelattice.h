
#ifndef INITIALIZELATTICE_H
#define INITIALIZELATTICE_H

void InitializeLattice ( int dim, double **SpinMatrix, double &E, double &M, int ordered );

//function for indexing periodic boundary conditions
inline int PeriodicBoundary( int i, int limit, int add ) {
    return ( i + limit + add ) % ( limit );
}

#endif // INITIALIZELATTICE_H

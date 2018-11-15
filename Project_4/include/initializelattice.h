
#ifndef INITIALIZELATTICE_H
#define INITIALIZELATTICE_H

// inline int PeriodicBoundary( int i, int limit, int add );
void InitializeLattice ( int dim, double **SpinMatrix, double &E, double &M, int ordered );

inline int PeriodicBoundary( int i, int limit, int add ) {
    return ( i + limit + add ) % ( limit );
}

#endif // INITIALIZELATTICE_H

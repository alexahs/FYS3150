
#ifndef METROPOLISSAMPLING_H
#define METROPOLISSAMPLING_H

void MetropolisSampling ( int dim, int MCcycles,int loopStart, int loopStop, double T, double *ExpectVal,
                          int ordered, double &Accept);

#endif // METROPOLISSAMPLING_H

#ifndef __MISAK_H__
#define __MISAK_H__

#include "TVector3.h"

void misak_init(); // Call to initialize misak's subroutines. 

double sigma_misak_df1(double Ebeam, TVector3 k, TVector3 p, bool isProton);  // Misak's implementation of CC1
double sigma_misak_st(double Ebeam, TVector3 k, TVector3 p, bool isProton);   // Stationary nucleon
double sigma_misak_free(double Ebeam, TVector3 k, TVector3 p, bool isProton); // Free (moving) nucleon
double sigma_misak_lc(double Ebeam, TVector3 k, TVector3 p, bool isProton);   // Full light-cone cross section

#endif

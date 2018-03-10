#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// Minimum recoil momenum
const double min_prec=0.35;

// For the binned acceptance correction
const int n_bins = 60;
const int n_slices = 56;
//const double samples_per_slice = 50000.; // Maybe this is a little high. Reduce to 10000?
const double samples_per_slice = 10000.; // Maybe this is a little high. Reduce to 10000?

// For the KDE acceptance correction
const double kde_width = 0.05;
const int n_pseudo=10;

// For determining the pp to p correction
const int n_pmiss_bins=28;

#endif

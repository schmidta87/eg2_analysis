#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// For the binned acceptance correction
const int n_bins = 60;
const int n_slices = 56;
const double samples_per_slice = 50000.;

// For the KDE acceptance correction
const double kde_width = 0.05;
const int n_pseudo=10;

#endif

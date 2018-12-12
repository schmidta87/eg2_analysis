#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// Minimum recoil momenum
const double min_prec=0.35;

// For the binned acceptance correction
const int n_bins = 60;
const int n_slices = 56;
//const double samples_per_slice = 50000.; // Maybe this is a little high. Reduce to 10000?
const double samples_per_slice = 10000.; 

// For the KDE acceptance correction
const double kde_width = 0.05;
const int n_pseudo=10;

// For determining the pp to p correction
const int n_pmiss_bins=28; // fine binning
const int n_pmiss_bins_coarse=10; // coarse binning
const double coarse_bin_edges[11]={0.3,0.375,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.875,1.0};


// Things for calculating cross sections, etc.
const double me = 0.000511;
const double mN = 0.93892;
const double mU=0.931;
const double GeVfm=0.1973;
const double alpha=0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;

#endif

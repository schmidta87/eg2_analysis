#include "fiducials.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "TVector3.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"

#include "constants.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Call fid_test using: \n"
	   << "\tfid_test /path/to/data/file /path/to/out/file\n\n";
      return -1;
    }
  
  ofstream outfile(argv[2]);
  TFile * infile = new TFile(argv[1]);
  TTree * t = (TTree*) infile->Get("T");
  Float_t Pe[3]; // Electron momenta
  Float_t phi_e;
  t->SetBranchAddress("Pe",Pe);
  t->SetBranchAddress("phi_e",&phi_e);

  int n_events = t->GetEntries();

  for (int i=0 ; i < n_events ; i++)
    {
      t->GetEntry(i);

      TVector3 p(Pe[0],Pe[1],Pe[2]);

      int fid = accept_electron(p)? 1: 0;
      
      outfile << p.Theta()*180./M_PI << " " << phi_e << " " << p.Mag() << " " << fid << "\n";
    }

  outfile.close();

  return 0;
}

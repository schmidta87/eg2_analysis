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
  const int maxN=10;
  Int_t nmb; // Number of protons in each event
  Float_t Pp[maxN][3]; // Proton momenta
  Float_t Rp[maxN][3]; // Proton vertices
  t->SetBranchAddress("nmb",&nmb);
  t->SetBranchAddress("Pp",Pp);
  t->SetBranchAddress("Rp",Rp);

  int n_events = t->GetEntries();

  for (int i=0 ; i < n_events ; i++)
    {
      t->GetEntry(i);

      for (int j=0 ; j<nmb ; j++)
	{
	  // Require coming from the solid target
	  double z = Rp[j][2];
	  if (fabs(z+22.25)>2.25)
	    continue;

	  TVector3 p(Pp[j][0],Pp[j][1],Pp[j][2]);
	  double sanPhiDeg = p.Phi() * 180./M_PI;
	  if (sanPhiDeg < -30.)
	    sanPhiDeg += 360.;

	  outfile << p.Mag() << " " << p.Theta()*180./M_PI << " " << sanPhiDeg << " " 
		  << accept_proton_simple(p) << " " << accept_proton(p) << "\n";
	}
    }

  outfile.close();

  return 0;
}

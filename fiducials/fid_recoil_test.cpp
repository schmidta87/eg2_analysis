#include "fiducials.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "TVector3.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Call fid_test using: \n"
	   << "\tfid_test /path/to/epp/file /path/to/out/file\n\n";
      return -1;
    }
  
  ofstream outfile(argv[2]);
  TFile * eppfile = new TFile(argv[1]);
  TTree * tepp = (TTree*) eppfile->Get("T");
  Float_t pPVec[2][3]; // Proton momenta
  Float_t pMissVec[2][3];
  Float_t vertices[2][3]; // position 3 vectors of proton vertices
  tepp->SetBranchAddress("Pmiss",pMissVec);
  tepp->SetBranchAddress("Rp",vertices);
  tepp->SetBranchAddress("Pp",pPVec);  
  int n_epp_events = tepp->GetEntries();

  for (int i=0 ; i < n_epp_events ; i++)
    {
      tepp->GetEntry(i);

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 prec(pPVec[1][0],pPVec[1][1],pPVec[1][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);

      // Make Or's cut on leading proton
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
        continue;

      // Make Or's cut on recoil proton
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > 0.35)))
        continue;

      int fid = accept_proton(prec)? 1: 0;
      
      outfile << prec.Theta()*180./M_PI << " " << prec.Phi()*180./M_PI << " " << prec.Mag() << " " << fid << "\n";
    }

  outfile.close();

  return 0;
}

#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorT.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "/path/to/2p/file /path/to/out/file\n\n";
      return -1;
    }

  TFile * outfile = new TFile(argv[2],"RECREATE");

  // Initialize some histograms
  TH2D * hist_epp_cm_lon = new TH2D("cm_lon","epp events;p_miss [GeV];p_cm_lon [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_inp = new TH2D("cm_inp","epp events;p_miss [GeV];p_cm_inp [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_oop = new TH2D("cm_oop","epp events;p_miss [GeV];p_cm_oop [GeV];Counts",7,0.3,1.,24,-1.2,1.2);

  // Read in epp file
  cerr << "Loading 2p data file " << argv[2] << " ...\n";
  Float_t pPVec[2][3]; // Proton momenta
  Float_t qVec[3]; // q
  Float_t pMissVec[2][3];
  Float_t vertices[2][3]; // position 3 vectors of proton vertices
  TFile *fepp = new TFile(argv[1]);
  TTree *tepp = (TTree*)fepp->Get("T");
  tepp->SetBranchAddress("q",qVec);
  tepp->SetBranchAddress("Pmiss",pMissVec);
  tepp->SetBranchAddress("Rp",vertices);
  tepp->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tepp->GetEntries() ; i++)
    {
      tepp->GetEvent(i);  

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 prec(pPVec[1][0],pPVec[1][1],pPVec[1][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);

      // Make Or's cut on leading proton
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
	continue;

      // Make Or's cut on recoil proton
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > 0.35)))
	continue;

      double pmiss_mag = pmiss.Mag();

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = lon.Cross(TVector3(0.,0.,1.)).Unit();
      TVector3 inp = lon.Cross(oop).Unit();

      TVector3 pcm = prec + pmiss;
      double pcm_lon = pcm.Dot(lon);
      double pcm_inp = pcm.Dot(inp);
      double pcm_oop = pcm.Dot(oop);

      hist_epp_cm_lon->Fill(pmiss_mag,pcm_lon);
      hist_epp_cm_inp->Fill(pmiss_mag,pcm_inp);
      hist_epp_cm_oop->Fill(pmiss_mag,pcm_oop);

    }
  fepp->Close();

  cerr << "Done reading input data.\n";

  // Write out histograms
  outfile->cd();
  hist_epp_cm_lon->Write();
  hist_epp_cm_inp->Write();
  hist_epp_cm_oop->Write();
  outfile->Close();

  // Clean up
  return 0;
}

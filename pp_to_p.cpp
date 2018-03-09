#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

#include "fiducials.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "/path/to/1p/file /path/to/2p/file /path/to/out/file\n\n";
      return -1;
    }

  TFile * outfile = new TFile(argv[3],"RECREATE");

  // Initialize some histograms
  double binEdges[11]={0.3,0.375,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.875,1.0};
  TH1D * hist_ep = new TH1D("ep","ep events;p_miss [GeV];Counts",10,binEdges);
  TH1D * hist_epp = new TH1D("epp","epp events;p_miss [GeV];Counts",10,binEdges);
  TH1D * hist_epp_nofid = new TH1D("epp_nofid","No fiducial cuts applied;p_miss [GeV];Counts",10,binEdges);
  TH1D * hist_ep_fine = new TH1D("ep_fine","ep events;p_miss [GeV];Counts",28,0.3,1.);
  TH1D * hist_epp_fine = new TH1D("epp_fine","epp events;p_miss [GeV];Counts",28,0.3,1.);
  TGraphAsymmErrors * pp_to_p = new TGraphAsymmErrors();
  pp_to_p->SetName("pp_to_p");
  pp_to_p->SetTitle("pp_to_p;p_miss [GeV];pp_to_p ratio");
  TGraphAsymmErrors * pp_to_p_nofid = new TGraphAsymmErrors();
  pp_to_p_nofid->SetName("pp_to_p_nofid");
  pp_to_p_nofid->SetTitle("No fiducial cuts;p_miss [GeV];pp_to_p ratio");

  // Read in epp file
  cerr << "Loading 2p data file " << argv[2] << " ...\n";
  Float_t pPVec[2][3]; // Proton momenta
  Float_t pMissVec[2][3];
  Float_t vertices[2][3]; // position 3 vectors of proton vertices
  TFile *fepp = new TFile(argv[2]);
  TTree *tepp = (TTree*)fepp->Get("T");
  tepp->SetBranchAddress("Pmiss",pMissVec);
  tepp->SetBranchAddress("Rp",vertices);
  tepp->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tepp->GetEntries() ; i++)
    {
      tepp->GetEvent(i);  

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 prec(pPVec[1][0],pPVec[1][1],pPVec[1][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      
      // Make Or's cut on leading proton
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
	continue;

      hist_ep->Fill(pmiss.Mag());
      hist_ep_fine->Fill(pmiss.Mag());

      // Make Or's cut on recoil proton
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > 0.35)))
	continue;

      hist_epp_nofid->Fill(pmiss.Mag());

      // Make fiducial cuts on the recoil proton
      if (!accept_proton(prec))
	continue;

      hist_epp->Fill(pmiss.Mag());
      hist_epp_fine->Fill(pmiss.Mag());
    }
  fepp->Close();

  // Read in ep file
  cerr << "Loading 1p data file " << argv[1] << " ...\n";
  TFile *fep = new TFile(argv[1]);
  TTree *tep = (TTree*)fep->Get("T");
  tep->SetBranchAddress("Pmiss",pMissVec);
  tep->SetBranchAddress("Rp",vertices);
  tep->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tep->GetEntries() ; i++)
    {
      tep->GetEvent(i);  

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 prec(pPVec[1][0],pPVec[1][1],pPVec[1][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      
      // Make Or's cut on leading proton
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
	continue;

      hist_ep->Fill(pmiss.Mag());
      hist_ep_fine->Fill(pmiss.Mag());

    }
  fep->Close();
  
  cerr << "Done reading input data.\n";
  
  pp_to_p->BayesDivide(hist_epp,hist_ep);
  pp_to_p_nofid->BayesDivide(hist_epp_nofid,hist_ep);
  
  // Write out histograms
  outfile->cd();
  hist_epp->Write();
  hist_ep->Write();
  hist_epp_fine->Write();
  hist_epp_nofid->Write();
  hist_ep_fine->Write();
  pp_to_p->Write();
  pp_to_p_nofid->Write();
  outfile->Close();

  // Clean up
  return 0;
}

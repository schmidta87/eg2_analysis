#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

#include "constants.h"
#include "AccMap.h"

using namespace std;

// Vectors which describe all of the useful info for the e'p events

const int n_recoils = 20;

int main(int argc, char **argv)
{
  if (argc != 9)
    {
      cerr << "Wrong number of arguments! Instead use\n"
	   << "\t /path/to/ep/file /path/to/epp/file /path/to/map/file [a1] [a2] [b1] [b2] [sigPerp] \n\n";
      return -1;
    }

  // Parameters
  const double a1=atof(argv[4]);
  const double a2=atof(argv[5]);
  const double b1=atof(argv[6]);
  const double b2=atof(argv[7]);
  const double sigPerp = atof(argv[8]);
  const double muPerp=0.;

  // Some set-up
  vector<TVector3> ep_pmiss_list;
  vector<TVector3> ep_q_list;
  vector<TVector3> ep_lon_list;
  vector<TVector3> ep_inp_list;
  vector<TVector3> ep_oop_list;
  TRandom3 myRand(0);

  // Read in e'p events
  cerr << "Loading 1p data file " << argv[1] << " ...\n";

  TFile * fep = new TFile(argv[1]);
  TTree * tep = (TTree*)fep->Get("T");
  Float_t pPVec[2][3]; // Proton momenta
  Float_t qVec[3]; // q
  Float_t pMissVec[2][3];
  Float_t vertices[2][3]; // position 3 vectors of proton vertices
  tep->SetBranchAddress("q",qVec);
  tep->SetBranchAddress("Pmiss",pMissVec);
  tep->SetBranchAddress("Rp",vertices);
  tep->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tep->GetEntries() ; i++)
    {
      tep->GetEvent(i);

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);

      // Make Or's cuts on vertex position
      if(!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
        continue;

      ep_pmiss_list.push_back(pmiss);
      ep_q_list.push_back(q);

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = oop.Cross(lon).Unit();
      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);
    }
  fep->Close();

  // Read in epp file
  cerr << "Loading 2p data file " << argv[2] << " ...\n";

  TFile *fepp = new TFile(argv[2]);
  TTree *tepp = (TTree*)fepp->Get("T");
  tepp->SetBranchAddress("q",qVec);
  tepp->SetBranchAddress("Pmiss",pMissVec);
  tepp->SetBranchAddress("Rp",vertices);
  tepp->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tepp->GetEntries() ; i++)
    {
      tepp->GetEvent(i);

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);

      // Make Or's cut on vertex position
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
        continue;

      ep_pmiss_list.push_back(pmiss);
      ep_q_list.push_back(q);
      double pmiss_mag = pmiss.Mag();

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = oop.Cross(lon).Unit();
      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);

    }
  fepp->Close();

  double n_ep_events = ep_pmiss_list.size();

  cerr << "Done reading input data.\n";
  cerr << "Read in " << n_ep_events << " ep events.\n";

  // Load the acceptance maps
  AccMap myMap(argv[3]);

  // Histograms for storing the acceptance data
  TH1D * h1Gen = new TH1D("gen","Generated;pmiss [GeV];Counts",n_pmiss_bins,0.3,1.);
  TH1D * h1Acc = new TH1D("acc","Accepted;pmiss [GeV];Counts",n_pmiss_bins,0.3,1.);
  h1Gen->Sumw2();
  h1Acc->Sumw2();
  TGraphAsymmErrors * corr = new TGraphAsymmErrors;

  // Generate a pseudo data sample for each e'p event
  // Loop over e'p
  for (int j=0 ; j< n_ep_events ; j++)
    {
      // Find pmiss, add to generated histogram
      double pmiss = ep_pmiss_list[j].Mag();
      h1Gen->Fill(pmiss,n_recoils);
      
      // Get the Longitudinal Gaussian parameters given pmiss
      double muLong = b1 * (pmiss - 0.6) + b2;
      double sigLong= a1 * (pmiss - 0.6) + a2;
      
      // Generate some p recoil vectors
      for (int k=0 ; k<n_recoils ; k++)
	{
	  double pcmLon = muLong + sigLong * myRand.Gaus();
	  double pcmInp = muPerp + sigPerp * myRand.Gaus();
	  double pcmOop = muPerp + sigPerp * myRand.Gaus();
	  TVector3 pcm = pcmLon * ep_lon_list[j] + pcmInp * ep_inp_list[j] + pcmOop * ep_oop_list[j];
	  TVector3 prec = pcm - ep_pmiss_list[j];
	  
	  // Test acceptance, add to accepted histogram
	  h1Acc->Fill(pmiss,myMap.recoil_accept(prec));
	}	  
    }
  
  corr->BayesDivide(h1Acc,h1Gen);
  // Loop over bins and fill to the 2d histogram
  for (int i=0 ; i< corr->GetN() ; i++)
    {
      double pmiss, correction;
      corr->GetPoint(i,pmiss,correction);

      cout << pmiss << " " << correction << " " << corr->GetErrorYlow(i) << " " << corr->GetErrorYhigh(i) << "\n";
    }
      
  return 0;
}

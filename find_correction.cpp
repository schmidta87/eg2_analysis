#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "constants.h"
#include "AccMap.h"

using namespace std;

// Vectors which describe all of the useful info for the e'p events

const int n_recoils = 20;

int main(int argc, char **argv)
{
  if (argc != 6)
    {
      cerr << "Wrong number of arguments! Instead use\n"
	   << "\t /path/to/ep/file /path/to/epp/file /path/to/map/file /path/to/mc/file /path/to/output/file \n\n";
      return -1;
    }

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

  // Create the input file and tree
  const double muPerp=0.;
  double a1, a2, b1, b2, sigPerp, logP, recoil_acc;
  TFile * infile = new TFile(argv[4]);
  TTree * intree = (TTree*) infile->Get("mcmc");
  intree->SetBranchAddress("a1",&a1);
  intree->SetBranchAddress("a2",&a2);
  intree->SetBranchAddress("b1",&b1);
  intree->SetBranchAddress("b2",&b2);
  intree->SetBranchAddress("sigPerp",&sigPerp);
  intree->SetBranchAddress("logposterior",&logP);
  intree->SetBranchAddress("recoil_acc",&recoil_acc);

  // Output file
  TFile * outfile = new TFile(argv[5],"RECREATE");
  TH2D * h2Corr = new TH2D("corr","Correction;pmiss [GeV];Acceptance [%];Counts",n_pmiss_bins,0.3,1.0,100,0.,100.);
  TH1D * h1Gen = new TH1D("gen","Generated;pmiss [GeV];Counts",n_pmiss_bins,0.3,1.);
  TH1D * h1Acc = new TH1D("acc","Accepted;pmiss [GeV];Counts",n_pmiss_bins,0.3,1.);
  // Loop over the tree
  for (int i=0 ; i<intree->GetEntries() ; i++)
    {
      if (i%100==0)
	cerr << "Working on event " << i << "\n";

      intree->GetEvent(i);

      // Clear the previous event's histograms
      h1Gen->Reset();
      h1Acc->Reset();

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

      // Loop over bins and fill to the 2d histogram
      for (int bin=1 ; bin <= n_pmiss_bins ; bin++)
	{
	  double pmiss = h1Gen->GetBinCenter(bin);
	  double gen = h1Gen->GetBinContent(bin);
	  double acc = h1Acc->GetBinContent(bin);

	  // sanitize in case no generated (maybe we can kill this)
	  double corr;
	  if (gen <= 0)
	    {
	      cerr << "MC sample " << i << ": no generated events for pmiss bin " << bin << " (pmiss = " << pmiss  << ")\n";
	      corr = 1.;
	    }
	  else
	    corr = acc/gen;

	  h2Corr->Fill(pmiss,corr*100); // Put in units of pct
	}
    }

  infile->Close();
      
  outfile->cd();
  h2Corr->Write();
  outfile->Close();

  return 0;
}

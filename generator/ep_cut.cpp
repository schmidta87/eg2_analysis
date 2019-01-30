#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include "AccMap.h"
#include "fiducials.h"
#include "Nuclear_Info.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tep_cut /path/to/gen/file /path/to/out/file\n\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  bool verbose = true;
  
  // Input Tree
  TTree * inTree = (TTree*)infile->Get("genT");
  Double_t gen_pe[3], gen_pLead[3], gen_pRec[3], gen_weight;
  Int_t lead_type, rec_type;
  inTree->SetBranchAddress("lead_type",&lead_type);
  inTree->SetBranchAddress("rec_type",&rec_type);
  inTree->SetBranchAddress("weight",&gen_weight);
  inTree->SetBranchAddress("pe",gen_pe);
  inTree->SetBranchAddress("pLead",gen_pLead);
  inTree->SetBranchAddress("pRec",gen_pRec);

  TFile * outfile = new TFile(argv[2],"RECREATE");

  // Output Tree
  TTree * outTree = new TTree("T","Cut Data Tree");
  Double_t pe[3], pLead[3], weight;
  outTree->Branch("weight",&weight,"weight/D");
  outTree->Branch("pe",&pe,"pe[3]/D");
  outTree->Branch("pLead",&pLead,"pLead[3]/D");

  // Loop over all events
  if (verbose)
    cerr << "Looping over 1p events...\n";
  const int nEvents = inTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %1000000==0 and verbose) 
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);

      // Require a leading proton
      if (lead_type != pCode)
	continue;

      // Create vectors for the particles
      TVector3 ve(gen_pe[0],gen_pe[1],gen_pe[2]);
      
      double gen_pe_Mag = ve.Mag();
      double gen_QSq = 2. * eg2beam * gen_pe_Mag * (1. - ve.CosTheta());
      double gen_nu = eg2beam - ve.Mag();
      double gen_xB = gen_QSq/(2.*mN*gen_nu);

      weight = gen_weight * 1.E33; // put it in nb to make it macroscopic

      if (weight <= 0.)
	continue;

      // Do leading proton cuts
      if (gen_xB < 1)
	continue;

      for (int i=0; i<3; i++)
	{
	  pe[i] = gen_pe[i];
	  pLead[i] = gen_pLead[i];
	}

      outTree->Fill();
    }

  outfile->cd();
  outTree->Write();
  outfile->Close();

  infile->Close();

  return 0;
}

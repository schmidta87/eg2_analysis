#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

#include "Nuclear_Info.h"
#include "fiducials.h"

using namespace std;

const double bin_edges[5] = {0.3,0.45,0.6,0.75,1.0};

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   pn_p_hists /path/to/input/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * fi = new TFile(argv[1]);
	TFile * fo = new TFile(argv[2],"RECREATE");

	// Create histograms        
	TH1D * hpn_pMiss = new TH1D("epn_Pmiss","epn;pMiss [GeV];Counts",35,0.3,1.0);
	hpn_pMiss->Sumw2();
	TH1D * hp_pMiss = new TH1D("ep_Pmiss","ep;pMiss [GeV];Counts",35,0.3,1.0);
	hp_pMiss->Sumw2();
	TH1D * hpn_pMiss_coarse = new TH1D("epn_Pmiss_coarse","epn;pMiss [GeV];Counts",4,bin_edges);
	hpn_pMiss_coarse->Sumw2();
	TH1D * hp_pMiss_coarse = new TH1D("ep_Pmiss_coarse","ep;pMiss [GeV];Counts",4,bin_edges);
	hp_pMiss_coarse->Sumw2();

	// pn2p graph
	TGraphAsymmErrors * pn_to_p = new TGraphAsymmErrors();
	pn_to_p->SetName("pn_to_p");
	pn_to_p->SetTitle("pn_to_p;p_miss [GeV]; pn_to_p ratio");
	TGraphAsymmErrors * pn_to_p_coarse = new TGraphAsymmErrors();
	pn_to_p_coarse->SetName("pn_to_p_coarse");
	pn_to_p_coarse->SetTitle("pn_to_p;p_miss [GeV]; pn_to_p ratio");
	
	// Loop over tree
	cerr << " Looping over tree...\n";
	TTree * ti = (TTree*)fi->Get("T");
	Double_t weightp, weightpn;
	Float_t Pmiss[2][3];
	Double_t pMiss;
	
	ti->SetBranchAddress("weightp",&weightp);
	ti->SetBranchAddress("weightpn",&weightpn);
	ti->SetBranchAddress("Pmiss",&Pmiss);
	
	for (int event =0 ; event < ti->GetEntries() ; event++)
	  {
	  
		ti->GetEvent(event);

		pMiss = sqrt(Pmiss[0][0]*Pmiss[0][0]+Pmiss[0][1]*Pmiss[0][1]+Pmiss[0][2]*Pmiss[0][2]);

		hp_pMiss->Fill(pMiss,weightp);
		hpn_pMiss->Fill(pMiss,weightpn);
		hp_pMiss_coarse->Fill(pMiss,weightp);
		hpn_pMiss_coarse->Fill(pMiss,weightpn);
	  }
	
	fi->Close();

	pn_to_p->Divide(hpn_pMiss,hp_pMiss,"cl=0.683 b(1,1) mode");
	pn_to_p_coarse->Divide(hpn_pMiss_coarse,hp_pMiss_coarse,"cl=0.683 b(1,1) mode");
	
	// Write out
	fo->cd();
	pn_to_p->Write();
	hpn_pMiss->Write();
	hp_pMiss->Write();
	pn_to_p_coarse->Write();
	hpn_pMiss_coarse->Write();
	hp_pMiss_coarse->Write();
	fo->Close();

	return 0;
}

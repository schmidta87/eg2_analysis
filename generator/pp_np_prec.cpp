#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TRatioPlot.h"

#include "Nuclear_Info.h"
#include "fiducials.h"

using namespace std;

int main(int argc, char ** argv)
{
	if (argc != 4)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   make_hists /path/to/pp/file /path/to/np/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * fpp = new TFile(argv[1]);
	TFile * fnp = new TFile(argv[2]);
	TFile * fo = new TFile(argv[3],"RECREATE");

	// Create histograms
        
	TH1D * hpp_pRec = new TH1D("epp_Pr","epp;pRec [GeV];Counts",26,0.35,1.0);
	hpp_pRec->Sumw2();
	TH1D * hnp_pRec = new TH1D("enp_Pr","enp;pRec [GeV];Counts",26,0.35,1.0);
	hnp_pRec->Sumw2();
	
	// Loop over pp tree
	cerr << " Looping over pp tree...\n";
	TTree * tpp = (TTree*)fpp->Get("T");
	Double_t weight;
	Double_t pRec;
        
	tpp->SetBranchAddress("weight",&weight);
	tpp->SetBranchAddress("pRec",&pRec);
        
	for (int event =0 ; event < tpp->GetEntries() ; event++)
	{
		tpp->GetEvent(event);

	        hpp_pRec->Fill(pRec,weight);
	}

	// Loop over np tree
	cerr << " Looping over np tree...\n";
	TTree * tnp = (TTree*)fnp->Get("T");

	tnp->SetBranchAddress("weight",&weight);
	tnp->SetBranchAddress("pRec",&pRec);
        
      	for (int event =0 ; event < tnp->GetEntries() ; event++)
	{
		tnp->GetEvent(event);

	        hnp_pRec->Fill(pRec,weight);
	}
	
	fpp->Close();
	fnp->Close();

	// pp2np hist
	//TRatioPlot * pp_to_np = new TRatioPlot(hpp_pRec,hnp_pRec);

	TH1F *pp_to_np = (TH1F*)hpp_pRec->Clone("pp_to_np");
	pp_to_np->SetName("pp_to_np");
	pp_to_np->SetTitle("pp_to_np;p_rec [GeV];pp_to_np ratio");
	pp_to_np->Divide(hnp_pRec);
	pp_to_np->Scale(0.2);
	
	// Write out
	fo->cd();
	pp_to_np->Write();
	hpp_pRec->Write();
	hnp_pRec->Write();
	fo->Close();

	return 0;
}

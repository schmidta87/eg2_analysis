#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "Nuclear_Info.h"
#include "fiducials.h"

using namespace std;

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   prec_hists /path/to/input/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * fi = new TFile(argv[1]);;
	TFile * fo = new TFile(argv[2],"RECREATE");

	// Create histograms        
	//TH1D * hpp_pRec = new TH1D("epp_Pr","epp;pRec [GeV];Counts",26,0.35,1.0);
	//hpp_pRec->Sumw2();
	//TH1D * hnp_pRec = new TH1D("enp_Pr","enp;pRec [GeV];Counts",26,0.35,1.0);
	//hnp_pRec->Sumw2();

	Double_t bins[3] = {0.35,0.6,1.0};
	TH1D * hpp_pRec_2 = new TH1D("epp_Pr_2","epp;pRec [GeV];Counts",2,bins);
	hpp_pRec_2->Sumw2();
	TH1D * hnp_pRec_2 = new TH1D("enp_Pr_2","enp;pRec [GeV];Counts",2,bins);
	hnp_pRec_2->Sumw2();

	//TH1D * hpp_pMiss = new TH1D("epp_Pmiss","epp;pMiss [GeV];Counts",26,0.35,1.0);
	//hpp_pMiss->Sumw2();
	//TH1D * hnp_pMiss = new TH1D("enp_Pmiss","enp;pMiss [GeV];Counts",26,0.35,1.0);
	//hnp_pMiss->Sumw2();

	
	// Loop over pp tree
	cerr << " Looping over tree...\n";
	TTree * ti = (TTree*)fi->Get("T");
	Double_t weight;
	Double_t pRec;
	Int_t nump;
	//Double_t pLead;
	//Float_t Pmiss[2][3];
	//Double_t pMiss;
	
	ti->SetBranchAddress("weight",&weight);
	ti->SetBranchAddress("pRec",&pRec);
	ti->SetBranchAddress("nump",&nump);
	//ti->SetBranchAddress("pLead",&pLead);
	//ti->SetBranchAddress("Pmiss",&Pmiss);
	
	for (int event =0 ; event < ti->GetEntries() ; event++)
	  {
	  
		ti->GetEvent(event);

		//pMiss = sqrt(Pmiss[0][0]*Pmiss[0][0]+Pmiss[0][1]*Pmiss[0][1]+Pmiss[0][2]*Pmiss[0][2]);

		if (nump == 2)
		  {
		    //hpp_pRec->Fill(pRec,weight);
		    hpp_pRec_2->Fill(pRec,weight);
		    
		    //hpp_pMiss->Fill(pMiss,weight);
		  }
		else
		  {
		    //hnp_pRec->Fill(pRec,weight);
		    hnp_pRec_2->Fill(pRec,weight);

		    //hnp_pMiss->Fill(pMiss,weight);
		  }
	  }
	
	fi->Close();

	// pp2np hist
	//TH1F *pp_to_np = (TH1F*)hpp_pRec->Clone("pp_to_np");
	//pp_to_np->SetName("pp_to_np");
	//pp_to_np->SetTitle("pp_to_np;p_rec [GeV];pp_to_np ratio");
	//pp_to_np->Divide(hnp_pRec);
	
	TH1F *pp_to_np_2 = (TH1F*)hpp_pRec_2->Clone("pp_to_np_2");
	pp_to_np_2->SetName("pp_to_np_2");
	pp_to_np_2->SetTitle("pp_to_np_2;p_rec [GeV];pp_to_np ratio");
	pp_to_np_2->Divide(hnp_pRec_2);	

	//TH1F *pp_to_np_miss = (TH1F*)hpp_pMiss->Clone("pp_to_np_miss");
	//pp_to_np_miss->SetName("pp_to_np_miss");
	//pp_to_np_miss->SetTitle("pp_to_np_miss;p_miss [GeV];pp_to_np ratio");
	//pp_to_np_miss->Divide(hnp_pMiss);
	
	//pp_to_np->SetMaximum(0.4);
	//pp_to_np->SetMinimum(0.0);
	pp_to_np_2->SetMaximum(0.4);
	pp_to_np_2->SetMinimum(0.0);
	//hpp_pRec->SetMinimum(0.0);
	//hnp_pRec->SetMinimum(0.0);
	//hpp_pRec_2->SetMinimum(0.0);
	//hnp_pRec_2->SetMinimum(0.0);
        
	// Write out
	fo->cd();
	//pp_to_np->Write();
	//hpp_pRec->Write();
	//hnp_pRec->Write();
	pp_to_np_2->Write();
	hpp_pRec_2->Write();
	hnp_pRec_2->Write();
	//pp_to_np_miss->Write();
	//hpp_pMiss->Write();
	//hnp_pMiss->Write();
	fo->Close();

	return 0;
}
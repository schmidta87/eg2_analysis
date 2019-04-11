#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "Nuclear_Info.h"
#include "fiducials.h"
#include "helpers.h"

using namespace std;

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   pp_np_likelihood /path/to/input/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * fi = new TFile(argv[1]);

	Double_t bins[3] = {0.35,0.6,1.0};
	TH1D * hpp_pRec_2 = new TH1D("epp_Pr_2","epp;pRec [GeV];Counts",2,bins);
	hpp_pRec_2->Sumw2();
	TH1D * hnp_pRec_2 = new TH1D("enp_Pr_2","enp;pRec [GeV];Counts",2,bins);
	hnp_pRec_2->Sumw2();
	
	// Loop over pp tree
	cerr << " Looping over tree...\n";
	TTree * ti = (TTree*)fi->Get("T");
	Double_t weight;
	Double_t pRec;
	Int_t nump;
	
	ti->SetBranchAddress("weight",&weight);
	ti->SetBranchAddress("pRec",&pRec);
	ti->SetBranchAddress("nump",&nump);

	for (int event =0 ; event < ti->GetEntries() ; event++)
	  {

		ti->GetEvent(event);

		if (nump == 2)
		  {
		    hpp_pRec_2->Fill(pRec,weight);		    
		  }
		else
		  {
		    hnp_pRec_2->Fill(pRec,weight);
		  }
	  }


	TH1F *pp_to_np_2 = (TH1F*)hpp_pRec_2->Clone("pp_to_np_2");
	pp_to_np_2->SetName("pp_to_np_2");
	pp_to_np_2->SetTitle("pp_to_np_2;p_rec [GeV];pp_to_np ratio");
	pp_to_np_2->Divide(hnp_pRec_2);	

	// Write out

	double likelihood = exp(-sq(pp_to_np_2->GetBinContent(1)*100-5.33)/(2*sq(0.74))
				-sq(pp_to_np_2->GetBinContent(2)*100-7.87)/(2*sq(1.82)));

	fi->Close();

	std::ofstream ofs (argv[2], std::ofstream::app);
	ofs << likelihood << "\n";
	ofs.close();

	return 0;
}

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TVectorT.h"

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
	TH1D * hpn_Pmr = new TH1D("epn_Pmr","epn;Theta_Pmr [deg];Counts",20,100.,180.);
	hpn_Pmr->Sumw2();
	TH1D * hpn_cPmr = new TH1D("epn_cPmr","epn;cos(Theta_Pmr);Counts",20,-1.,0.);
	hpn_cPmr->Sumw2();
	TH1D * hp_Emiss_fine = new TH1D("ep_Emiss_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	hp_Emiss_fine->Sumw2();
	TH1D * hpn_Emiss_fine = new TH1D("epn_Emiss_fine","epn;Emiss [GeV];Counts",160,-0.2,0.6);
	hpn_Emiss_fine->Sumw2();

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
	Float_t Q2, Xb, Pmiss[2][3], Pmiss_size[2], Pp_size[2];
	Double_t pMiss;
	Float_t Pp[2][3], Pe[3], q[3];
	
	ti->SetBranchAddress("weightp",&weightp);
	ti->SetBranchAddress("weightpn",&weightpn);
	ti->SetBranchAddress("Xb",&Xb);
	ti->SetBranchAddress("Q2",&Q2);
	ti->SetBranchAddress("Pmiss",&Pmiss);
	ti->SetBranchAddress("Pmiss_size",&Pmiss_size);
	ti->SetBranchAddress("Pp_size",&Pp_size);
	ti->SetBranchAddress("Pp",Pp);
	ti->SetBranchAddress("Pe",&Pe);
	ti->SetBranchAddress("q",&q);
	
	for (int event =0 ; event < ti->GetEntries() ; event++)
	  {
	  
		ti->GetEvent(event);
		
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);  
		TVector3 vq(q[0],q[1],q[2]);
		TVector3 vmiss = vlead - vq;
		TVector3 vcm = vmiss + vrec; 

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);
		pMiss = vmiss.Mag();

		hp_pMiss->Fill(pMiss,weightp);
		hpn_pMiss->Fill(pMiss,weightpn);
		hp_pMiss_coarse->Fill(pMiss,weightp);
		hpn_pMiss_coarse->Fill(pMiss,weightpn);
		hpn_Pmr->Fill(vmiss.Angle(vrec)*180./M_PI,weightpn);
		hpn_cPmr->Fill(cos(vmiss.Angle(vrec)),weightpn);
		hp_Emiss_fine->Fill(Emiss,weightp);
		hpn_Emiss_fine->Fill(Emiss,weightpn);
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
	hpn_Pmr->Write();
	hpn_cPmr->Write();
	hp_Emiss_fine->Write();
	hpn_Emiss_fine->Write();


	fo->Close();

	return 0;
}

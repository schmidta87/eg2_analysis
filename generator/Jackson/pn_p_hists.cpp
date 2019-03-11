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

const double pmiss_lo=0.45;
const double pmiss_md=0.6;
const double pmiss_hi=0.75;

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

	vector<TH1*> h_list;

	// Create histograms        
	TH1D * hpn_pMiss = new TH1D("epn_Pmiss","epn;pMiss [GeV];Counts",35,0.3,1.0);
	h_list.push_back(hpn_pMiss);
	TH1D * hp_pMiss = new TH1D("ep_Pmiss","ep;pMiss [GeV];Counts",35,0.3,1.0);
	h_list.push_back(hp_pMiss);
	TH1D * hpn_pMiss_tot = new TH1D("epn_Pmiss_tot","epn;pMiss [GeV];Counts",1,0.3,1.0);
	h_list.push_back(hpn_pMiss_tot);
	TH1D * hp_pMiss_tot = new TH1D("ep_Pmiss_tot","ep;pMiss [GeV];Counts",1,0.3,1.0);
	h_list.push_back(hp_pMiss_tot);
	TH1D * hpn_pMiss_coarse = new TH1D("epn_Pmiss_coarse","epn;pMiss [GeV];Counts",4,bin_edges);
	h_list.push_back(hpn_pMiss_coarse);
	TH1D * hp_pMiss_coarse = new TH1D("ep_Pmiss_coarse","ep;pMiss [GeV];Counts",4,bin_edges);
	h_list.push_back(hp_pMiss_coarse);
	TH1D * hpn_Pmr = new TH1D("epn_Pmr","epn;Theta_Pmr [deg];Counts",20,100.,180.);
	h_list.push_back(hpn_Pmr);
	TH1D * hpn_cPmr = new TH1D("epn_cPmr","epn;cos(Theta_Pmr);Counts",20,-1.,0.);
	h_list.push_back(hpn_cPmr);
	TH1D * hp_Emiss_fine = new TH1D("ep_Emiss_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h_list.push_back(hp_Emiss_fine);
	TH1D * hpn_Emiss_fine = new TH1D("epn_Emiss_fine","epn;Emiss [GeV];Counts",160,-0.2,0.6);
	h_list.push_back(hpn_Emiss_fine);
	TH1D * hpn_Emiss = new TH1D("epn_Emiss","epn;Emiss [GeV];Counts",20,-0.1,0.5);
	h_list.push_back(hpn_Emiss);
	TH1D * hpn_mom2 = new TH1D("epn_mom2","epp;Recoil Mom [GeV/c];Counts",17,0.35,1.2);
	h_list.push_back(hpn_mom2);

	TH1D * hpn_QSq_split[4];
	TH1D * hpn_mMiss_split[4];
	TH1D * hpn_xB_split[4];
	TH1D * hpn_q_split[4];
	TH1D * hpn_Emiss_split[4];

	for (int i=0; i<4; i++)
	  {
	    char temp[100];
	    
	    sprintf(temp,"epn_QSq_%d",i);
	    hpn_QSq_split[i] = new TH1D(temp,"ep;QSq [GeV^2];Counts",20,1.,3.5);
	    h_list.push_back(hpn_QSq_split[i]);
	    
	    sprintf(temp,"epn_mMiss_%d",i);
	    hpn_mMiss_split[i] = new TH1D(temp,"ep;mMiss [GeV];Counts",20,0.6,1.2);
	    h_list.push_back(hpn_mMiss_split[i]);
	    
	    sprintf(temp,"epn_xB_%d",i);
	    hpn_xB_split[i] = new TH1D(temp,"ep;xB;Counts",20,1.,2.5);
	    h_list.push_back(hpn_xB_split[i]);
	    
	    sprintf(temp,"epn_q_%d",i);
	    hpn_q_split[i] = new TH1D(temp,"ep;q [GeV];Counts",20,1.,2.5);
	    h_list.push_back(hpn_q_split[i]);
	    
	    sprintf(temp,"epn_Emiss_%d",i);
	    hpn_Emiss_split[i] = new TH1D(temp,"ep;Emiss [GeV];Counts",20,-0.1,0.5);
	    h_list.push_back(hpn_Emiss_split[i]);
	  }

	for (int i=0; i<h_list.size(); i++)
	  h_list[i]->Sumw2();

	// pn2p graph
	TGraphAsymmErrors * pn_to_p = new TGraphAsymmErrors();
	pn_to_p->SetName("pn_to_p");
	pn_to_p->SetTitle("pn_to_p;p_miss [GeV]; pn_to_p ratio");
	TGraphAsymmErrors * pn_to_p_coarse = new TGraphAsymmErrors();
	pn_to_p_coarse->SetName("pn_to_p_coarse");
	pn_to_p_coarse->SetTitle("pn_to_p;p_miss [GeV]; pn_to_p ratio");
	TGraphAsymmErrors * pn_to_p_tot = new TGraphAsymmErrors();
	pn_to_p_tot->SetName("pn_to_p_tot");
	pn_to_p_tot->SetTitle("pn_to_p;p_miss [GeV]; pn_to_p ratio");
	
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

		pMiss = vmiss.Mag();
		double omega = Q2/(2.*mN*Xb);
		double ELead = sqrt(mN*mN+vlead.Mag2());
		double m_miss = sqrt((omega+2*mN-ELead)*(omega+2*mN-ELead)-vmiss.Mag2());
		double Emiss = -m_12C + mN + sqrt((omega + m_12C - ELead)*(omega + m_12C - ELead) - (Pmiss_size[0]*Pmiss_size[0]));

		hp_pMiss->Fill(pMiss,weightp);
		hpn_pMiss->Fill(pMiss,weightpn);
		hp_pMiss_coarse->Fill(pMiss,weightp);
		hpn_pMiss_coarse->Fill(pMiss,weightpn);
		hp_pMiss_tot->Fill(pMiss,weightp);
		hpn_pMiss_tot->Fill(pMiss,weightpn);
		hpn_Pmr->Fill(vmiss.Angle(vrec)*180./M_PI,weightpn);
		hpn_cPmr->Fill(cos(vmiss.Angle(vrec)),weightpn);
		hp_Emiss_fine->Fill(Emiss,weightp);
		hpn_Emiss_fine->Fill(Emiss,weightpn);
		hpn_Emiss->Fill(Emiss,weightpn);


		int Pmiss_region;
		if (Pmiss_size[0] < pmiss_lo) Pmiss_region = 0;
		else if (Pmiss_size[0] < pmiss_md) Pmiss_region = 1;
		else if (Pmiss_size[0] < pmiss_hi) Pmiss_region = 2;
		else Pmiss_region = 3;
		
		hpn_QSq_split[Pmiss_region]->Fill(Q2,weightpn);
		hpn_xB_split[Pmiss_region]->Fill(Xb,weightpn);
		hpn_q_split[Pmiss_region]->Fill(vq.Mag(),weightpn);
		hpn_mMiss_split[Pmiss_region]->Fill(m_miss,weightpn);
		hpn_Emiss_split[Pmiss_region]->Fill(Emiss,weightpn);

		hpn_mom2->Fill(vrec.Mag(),weightpn);
		
	  }
	
	fi->Close();

	pn_to_p->Divide(hpn_pMiss,hp_pMiss,"cl=0.683 b(1,1) mode");
	pn_to_p_coarse->Divide(hpn_pMiss_coarse,hp_pMiss_coarse,"cl=0.683 b(1,1) mode");
	pn_to_p_tot->Divide(hpn_pMiss_tot,hp_pMiss_tot,"cl=0.683 b(1,1) mode");
	
	// Write out
	fo->cd();

	pn_to_p->Write();
	pn_to_p_coarse->Write();
	pn_to_p_tot->Write();

	const double data_epn = 138.;
	const double norm = data_epn/hpn_pMiss->Integral();

	for (int i=0; i<h_list.size(); i++)
	  {
	    h_list[i]->Scale(norm);
	    h_list[i]->Write();
	  }

	fo->Close();

	return 0;
}

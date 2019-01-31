#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TVectorT.h"

#include "Nuclear_Info.h"
#include "fiducials.h"
#include "helpers.h"

using namespace std;

const double pmiss_lo=0.45;
const double pmiss_md=0.6;
const double pmiss_hi=0.75;

int main(int argc, char ** argv)
{
	if (argc != 4)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   make_hists /path/to/1p/file /path/to/2p/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * f1p = new TFile(argv[1]);
	TFile * f2p = new TFile(argv[2]);
	TFile * fo = new TFile(argv[3],"RECREATE");

	// Let's create a vector of all the histogram pointers so we can loop over them, save hassles
	vector<TH1*> h1p_list;
	vector<TH1*> h2p_list;

	// Create histograms
	TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",40,1.,5.);
	h1p_list.push_back(h1p_QSq);
	TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",26,1.2,2.5);
	h1p_list.push_back(h1p_xB );
	TH1D * h1p_Pm =  new TH1D("ep_Pm" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm );
	TH1D * h1p_Pm_coarse =  new TH1D("ep_Pm_coarse" ,"ep;pMiss [GeV];Counts",10,coarse_bin_edges);
	h1p_list.push_back(h1p_Pm_coarse);
	TH1D * h1p_Pmq = new TH1D("ep_Pmq","ep;Theta_Pmq [deg];Counts",40,100.,180.);
	h1p_list.push_back(h1p_Pmq);
	TH1D * h1p_cPmq = new TH1D("ep_cPmq","ep;cos(Theta_Pmq);Counts",40,-1.,0.);
	h1p_list.push_back(h1p_cPmq);
	TH1D * h1p_phi1 = new TH1D("ep_phi1","ep;Phi_1 [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phi1);
	TH1D * h1p_phie = new TH1D("ep_phie","ep;Phi_e [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phie);
	TH1D * h1p_theta1 = new TH1D("ep_theta1","ep;Theta_1 [deg];Counts",60,10.,130.);
	h1p_list.push_back(h1p_theta1);
	TH1D * h1p_thetae = new TH1D("ep_thetae","ep;Theta_e [deg];Counts",60,10.,40.);
	h1p_list.push_back(h1p_thetae);
	TH1D * h1p_mome = new TH1D("ep_mome","ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h1p_list.push_back(h1p_mome);
	TH1D * h1p_mom1 = new TH1D("ep_mom1","ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h1p_list.push_back(h1p_mom1);
	TH1D * h1p_Emiss = new TH1D("ep_Emiss","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss);
	TH1D * h1p_Emiss_lo = new TH1D("ep_Emiss_lo","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_lo);
	TH1D * h1p_Emiss_md = new TH1D("ep_Emiss_md","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_md);
	TH1D * h1p_Emiss_hi = new TH1D("ep_Emiss_hi","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi);
	TH1D * h1p_Emiss_hi1 = new TH1D("ep_Emiss_hi1","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi1);
	TH1D * h1p_Emiss_hi2 = new TH1D("ep_Emiss_hi2","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi2);
	TH1D * h1p_Tmiss_lo = new TH1D("ep_Tmiss_lo","ep;Tmiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Tmiss_lo);
	TH1D * h1p_Tmiss_md = new TH1D("ep_Tmiss_md","ep;Tmiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Tmiss_md);
	TH1D * h1p_Tmiss_hi = new TH1D("ep_Tmiss_hi","ep;Tmiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Tmiss_hi);
	TH1D * h1p_Tmiss_hi1 = new TH1D("ep_Tmiss_hi1","ep;Tmiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Tmiss_hi1);
	TH1D * h1p_Tmiss_hi2 = new TH1D("ep_Tmiss_hi2","ep;Tmiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Tmiss_hi2);
	TH1D * h1p_Emiss_fine = new TH1D("ep_Emiss_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_fine);
	TH1D * h1p_Emiss_lo_fine = new TH1D("ep_Emiss_lo_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_lo_fine);
	TH1D * h1p_Emiss_md_fine = new TH1D("ep_Emiss_md_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_md_fine);
	TH1D * h1p_Emiss_hi_fine = new TH1D("ep_Emiss_hi_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi_fine);
	TH1D * h1p_Emiss_hi1_fine = new TH1D("ep_Emiss_hi1_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi1_fine);
	TH1D * h1p_Emiss_hi2_fine = new TH1D("ep_Emiss_hi2_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi2_fine);
	TH2D * h1p_pmiss_Emiss = new TH2D("ep_pmiss_Emiss","ep;pmiss [GeV];Emiss [GeV];Counts",28,0.3,1.0,20,-0.2,0.6);
	h1p_list.push_back(h1p_pmiss_Emiss);
	TH2D * h1p_pmiss_E1 = new TH2D("ep_pmiss_E1","ep;pmiss [GeV];E1 [GeV];Counts",28,0.3,1.0,25,0.5,1.0);
	h1p_list.push_back(h1p_pmiss_E1);
	TH1D * h2p_QSq = new TH1D("epp_QSq","epp;QSq [GeV^2];Counts",40,1.,5.);
	h2p_list.push_back(h2p_QSq);
	TH1D * h2p_xB =  new TH1D("epp_xB" ,"epp;xB;Counts",26,1.2,2.5);
	h2p_list.push_back(h2p_xB );
	TH1D * h2p_Pm =  new TH1D("epp_Pm" ,"epp;pMiss [GeV];Counts",28,0.3,1.0);
	h2p_list.push_back(h2p_Pm );
	TH1D * h2p_Pm_coarse =  new TH1D("epp_Pm_coarse" ,"epp;pMiss [GeV];Counts",10,coarse_bin_edges);
	h2p_list.push_back(h2p_Pm_coarse);
	TH1D * h2p_Pmq = new TH1D("epp_Pmq","epp;Theta_Pmq [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmq);
	TH1D * h2p_cPmq = new TH1D("epp_cPmq","epp;cos(Theta_Pmq);Counts",20,-1.,0.);
	h2p_list.push_back(h2p_cPmq);
	TH1D * h2p_Pmr = new TH1D("epp_Pmr","epp;Theta_Pmr [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmr);
	TH1D * h2p_cPmr = new TH1D("epp_cPmr","epp;cos(Theta_Pmr);Counts",20,-1.,0.);
	h2p_list.push_back(h2p_cPmr);
	TH1D * h2p_phi1 = new TH1D("epp_phi1","epp;Phi_1 [deg];Counts",60,-30.,330.);
	h2p_list.push_back(h2p_phi1);
	TH1D * h2p_phi2 = new TH1D("epp_phi2","epp;Phi_2 [deg];Counts",60,-30.,330.);
	h2p_list.push_back(h2p_phi2);
	TH1D * h2p_phie = new TH1D("epp_phie","epp;Phi_e [deg];Counts",60,-30.,330.);
	h2p_list.push_back(h2p_phie);
	TH1D * h2p_thetae = new TH1D("epp_thetae","epp;Theta_e [deg];Counts",30,10.,40.);
	h2p_list.push_back(h2p_thetae);
	TH1D * h2p_mome = new TH1D("epp_mome","epp;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h2p_list.push_back(h2p_mome);
	TH1D * h2p_theta1 = new TH1D("epp_theta1","epp;Theta_1 [deg];Counts",30,10.,130.);
	h2p_list.push_back(h2p_theta1);
	TH1D * h2p_theta2 = new TH1D("epp_theta2","epp;Theta_2 [deg];Counts",30,10.,130.);
	h2p_list.push_back(h2p_theta2);
	TH1D * h2p_mom2 = new TH1D("epp_mom2","epp;Recoil Mom [GeV/c];Counts",17,0.35,1.2);
	h2p_list.push_back(h2p_mom2);
	TH1D * h2p_mom1 = new TH1D("epp_mom1","epp;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h2p_list.push_back(h2p_mom1);
	TH1D * h2p_Emiss = new TH1D("epp_Emiss","epp;Emiss [GeV];Counts",40,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss);
	TH1D * h2p_Emiss_lo = new TH1D("epp_Emiss_lo","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_lo);
	TH1D * h2p_Emiss_md = new TH1D("epp_Emiss_md","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_md);
	TH1D * h2p_Emiss_hi = new TH1D("epp_Emiss_hi","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi);
	TH1D * h2p_Emiss_hi1 = new TH1D("epp_Emiss_hi1","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi1);
	TH1D * h2p_Emiss_hi2 = new TH1D("epp_Emiss_hi2","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi2);
	TH1D * h2p_Emiss_fine = new TH1D("epp_Emiss_fine","epp;Emiss [GeV];Counts",160,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_fine);
	TH1D * h2p_Emiss_lo_fine = new TH1D("epp_Emiss_lo_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_lo_fine);
	TH1D * h2p_Emiss_md_fine = new TH1D("epp_Emiss_md_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_md_fine);
	TH1D * h2p_Emiss_hi_fine = new TH1D("epp_Emiss_hi_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi_fine);
	TH1D * h2p_Emiss_hi1_fine = new TH1D("epp_Emiss_hi1_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi1_fine);
	TH1D * h2p_Emiss_hi2_fine = new TH1D("epp_Emiss_hi2_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi2_fine);
	TH2D * h2p_pmiss_E1 = new TH2D("epp_pmiss_E1","epp;pmiss [GeV];E1 [GeV];Counts",28,0.3,1.0,25,0.5,1.0);
	h2p_list.push_back(h2p_pmiss_E1);
	TH2D * h2p_pmiss_appEstar = new TH2D("epp_pmiss_appEstar","epp;pmiss [GeV];Apparent Estar [GeV];Counts",28,0.3,1.0,20,-0.2,0.8);
	h2p_list.push_back(h2p_pmiss_appEstar);

	TH2D * pp_to_p_2d = new TH2D("pp_to_p_2d","2d ratio;pmiss [GeV];E1 [GeV];pp/p",28,0.3,1.0,20,0.5,0.9);

	TH1D * h1p_thetae_bySec[6];
	TH1D * h1p_theta1_bySec[6];
	TH1D * h2p_theta1_bySec[6];
	TH1D * h2p_theta2_bySec[6];
	for (int i=0 ; i<6 ; i++)
	{
		char temp[100];

		sprintf(temp,"ep_theta1_%d",i);
		h1p_theta1_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h1p_theta1_bySec[i]);

		sprintf(temp,"ep_thetae_%d",i);
		h1p_thetae_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,40.);
		h1p_list.push_back(h1p_thetae_bySec[i]);

		sprintf(temp,"epp_theta1_%d",i);
		h2p_theta1_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h2p_list.push_back(h2p_theta1_bySec[i]);

		sprintf(temp,"epp_theta2_%d",i);
		h2p_theta2_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h2p_list.push_back(h2p_theta2_bySec[i]);
	}

	// Now that all histograms have been defined, set them to Sumw2
	for (int i=0 ; i<h1p_list.size() ; i++)
	  h1p_list[i]->Sumw2();
	for (int i=0 ; i<h2p_list.size() ; i++)
	  h2p_list[i]->Sumw2();

	// pp2p graphs
	TGraphAsymmErrors * pp_to_p = new TGraphAsymmErrors();
	pp_to_p->SetName("pp_to_p");
	pp_to_p->SetTitle("pp_to_p;p_miss [GeV];pp_to_p ratio");
	TGraphAsymmErrors * pp_to_p_coarse = new TGraphAsymmErrors();
	pp_to_p_coarse->SetName("pp_to_p_coarse");
	pp_to_p_coarse->SetTitle("pp_to_p;p_miss [GeV];pp_to_p ratio");

	// Loop over 1p tree
	cerr << " Looping over 1p tree...\n";
	TTree * t1p = (TTree*)f1p->Get("T");
	Float_t Xb, Q2, Pmiss_size[2], Pp[2][3], Rp[2][3], Pp_size[2], Pmiss_q_angle[2], Pe[3], q[3];
	Double_t weight = 1.;
	t1p->SetBranchAddress("Pmiss_q_angle",Pmiss_q_angle);
	t1p->SetBranchAddress("Xb",&Xb);
	t1p->SetBranchAddress("Q2",&Q2);
	t1p->SetBranchAddress("Pmiss_size",Pmiss_size);
	t1p->SetBranchAddress("Pp_size",Pp_size);
	t1p->SetBranchAddress("Rp",Rp);
	t1p->SetBranchAddress("Pp",Pp);
	t1p->SetBranchAddress("Pe",Pe);
	// See if there is a weight branch
	TBranch * weight_branch = t1p->GetBranch("weight");
	if (weight_branch)
	{
		t1p->SetBranchAddress("weight",&weight);
	}

	for (int event =0 ; event < t1p->GetEntries() ; event++)
	{
		t1p->GetEvent(event);

		// Do necessary cuts
		if (fabs(Rp[0][2]+22.25)>2.25)
		  continue;
		if (Pp_size[0]>2.4)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;

		// Apply sector cuts
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = phie_deg/60.;
		if ((sec_e != 1) and (sec_e != 4) and (sec_e != 5))
		  continue;

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);
		double omega = Q2/(2.*mN*Xb);

		// Let's make a sanitized phi and sector
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome->Fill(ve.Mag(),weight);

		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = phi1_deg/60.;
		double theta1_deg = vp.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		
		// Let's figure out missing energy! 
		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);
		double Tmiss = sqrt(Pmiss_size[0]*Pmiss_size[0]+mN*mN)+m_10B-sqrt(Pmiss_size[0]*Pmiss_size[0]+m_11B*m_11B);
		h1p_Emiss->Fill(Emiss,weight);
		h1p_Emiss_fine->Fill(Emiss,weight);
		h1p_pmiss_Emiss->Fill(Pmiss_size[0],Emiss,weight);
		if (Pmiss_size[0] < pmiss_lo)
		  {
		    h1p_Emiss_lo->Fill(Emiss,weight);
		    h1p_Tmiss_lo->Fill(Tmiss,weight);
		    h1p_Emiss_lo_fine->Fill(Emiss,weight);
		  }
		else if (Pmiss_size[0] < pmiss_md)
		  {
		    h1p_Emiss_md->Fill(Emiss,weight);
		    h1p_Tmiss_md->Fill(Tmiss,weight);
		    h1p_Emiss_md_fine->Fill(Emiss,weight);
		  }
		else if (Pmiss_size[0] < 1.)
		  {
		    h1p_Emiss_hi->Fill(Emiss,weight);
		    h1p_Tmiss_hi->Fill(Tmiss,weight);
		    h1p_Emiss_hi_fine->Fill(Emiss,weight);
		    if (Pmiss_size[0] < pmiss_hi)
		      {
			h1p_Emiss_hi1->Fill(Emiss,weight);
			h1p_Tmiss_hi1->Fill(Tmiss,weight);
			h1p_Emiss_hi1_fine->Fill(Emiss,weight);
		      }
		    else
		      {
			h1p_Emiss_hi2->Fill(Emiss,weight);
			h1p_Tmiss_hi2->Fill(Tmiss,weight);
			h1p_Emiss_hi2_fine->Fill(Emiss,weight);
		      }
		  }

		h1p_pmiss_E1->Fill(Pmiss_size[0],sqrt(vp.Mag2() + mN*mN) - omega,weight);
	}

	// Loop over 2p tree
	cerr << " Looping over 2p tree...\n";
	TTree * t2p = (TTree*)f2p->Get("T");
	t2p->SetBranchAddress("Pmiss_q_angle",Pmiss_q_angle);
	t2p->SetBranchAddress("Xb",&Xb);
	t2p->SetBranchAddress("Q2",&Q2);
	t2p->SetBranchAddress("Pmiss_size",Pmiss_size);
	t2p->SetBranchAddress("Pp_size",Pp_size);
	t2p->SetBranchAddress("Rp",Rp);
	t2p->SetBranchAddress("Pp",Pp);
	t2p->SetBranchAddress("Pe",Pe);
	t2p->SetBranchAddress("q",q);
	// See if there is a weight branch
	weight=1.;
	weight_branch = t2p->GetBranch("weight");
	if (weight_branch)
	{
		t2p->SetBranchAddress("weight",&weight);
	}
	for (int event =0 ; event < t2p->GetEntries() ; event++)
	{
		t2p->GetEvent(event);

		// Do necessary cuts
		if (fabs(Rp[0][2]+22.25)>2.25)
			continue;
		if (Pp_size[0]>2.4)
			continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vlead))
		  continue;

		// Apply sector cuts
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = phie_deg/60.;
		if ((sec_e != 1) and (sec_e != 4) and (sec_e != 5))
		  continue;

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);
		double omega = Q2/(2.*mN*Xb);

		// Let's make a sanitized phi and sector
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome->Fill(ve.Mag(),weight);

		double phi1_deg = vlead.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = phi1_deg/60.;
		double theta1_deg = vlead.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1->Fill(Pp_size[0],weight);

		// Let's figure out missing energy! 
		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);
		double Tmiss = sqrt(Pmiss_size[0]*Pmiss_size[0]+mN*mN)+m_10B-sqrt(Pmiss_size[0]*Pmiss_size[0]+m_11B*m_11B);
		h1p_Emiss->Fill(Emiss,weight);
		h1p_Emiss_fine->Fill(Emiss,weight);
		h1p_pmiss_Emiss->Fill(Pmiss_size[0],Emiss,weight);
		h1p_pmiss_E1->Fill(Pmiss_size[0],sqrt(vlead.Mag2() + mN*mN) - omega,weight);
		if (Pmiss_size[0] < pmiss_lo)
		  {
		    h1p_Emiss_lo->Fill(Emiss,weight);
		    h1p_Tmiss_lo->Fill(Tmiss,weight);
		    h1p_Emiss_lo_fine->Fill(Emiss,weight);
		  }
		else if (Pmiss_size[0] < pmiss_md)
		  {
		    h1p_Emiss_md->Fill(Emiss,weight);
		    h1p_Tmiss_md->Fill(Tmiss,weight);
		    h1p_Emiss_md_fine->Fill(Emiss,weight);
		  }
		else if (Pmiss_size[0] < 1.)
		  {
		    h1p_Emiss_hi->Fill(Emiss,weight);
		    h1p_Tmiss_hi->Fill(Tmiss,weight);
		    h1p_Emiss_hi_fine->Fill(Emiss,weight);
		    if (Pmiss_size[0] < pmiss_hi)
		      {
			h1p_Emiss_hi1->Fill(Emiss,weight);
			h1p_Tmiss_hi1->Fill(Tmiss,weight);
			h1p_Emiss_hi1_fine->Fill(Emiss,weight);
		      }
		    else
		      {
			h1p_Emiss_hi2->Fill(Emiss,weight);
			h1p_Tmiss_hi2->Fill(Tmiss,weight);
			h1p_Emiss_hi2_fine->Fill(Emiss,weight);
		      }
		  }

		// Make a check on the recoils
		if (fabs(Rp[1][2]+22.25)>2.25)
		  continue;
		if (Pp_size[1] < 0.35)
		  continue;
		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);      
		if (!accept_proton(vrec))
		  continue;

		double Erec = sqrt(vrec.Mag2() + mN*mN);
		TVector3 vq(q[0],q[1],q[2]);
		TVector3 vmiss = vlead - vq;
		TVector3 vcm = vmiss + vrec;

		h2p_QSq->Fill(Q2,weight);
		h2p_xB ->Fill(Xb,weight);
		h2p_Pm ->Fill(Pmiss_size[0],weight);
		h2p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h2p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h2p_Pmr->Fill(vmiss.Angle(vrec)*180./M_PI,weight);
		h2p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);
		h2p_cPmr->Fill(cos(vmiss.Angle(vrec)),weight);

		h2p_phie->Fill(phie_deg,weight);
		h2p_thetae->Fill(thetae_deg,weight);
		h2p_mome->Fill(ve.Mag(),weight);

		h2p_phi1->Fill(phi1_deg,weight);
		h2p_theta1->Fill(theta1_deg,weight);
		h2p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h2p_mom1->Fill(Pp_size[0],weight);

		h2p_Emiss->Fill(Emiss,weight);
		h2p_Emiss_fine->Fill(Emiss,weight);
		h2p_pmiss_E1->Fill(Pmiss_size[0],sqrt(vlead.Mag2() + mN*mN) - omega,weight);

		if (Pmiss_size[0] < pmiss_lo)
		  {
		    h2p_Emiss_lo->Fill(Emiss,weight);
		    h2p_Emiss_lo_fine->Fill(Emiss,weight);
		  }
		else if (Pmiss_size[0] < pmiss_md)
		  {
		    h2p_Emiss_md->Fill(Emiss,weight);
		    h2p_Emiss_md_fine->Fill(Emiss,weight);
		  }
		else if (Pmiss_size[0] < 1.)
		  {
		    h2p_Emiss_hi->Fill(Emiss,weight);
		    h2p_Emiss_hi_fine->Fill(Emiss,weight);
		    if (Pmiss_size[0] < pmiss_hi)
		      {
			h2p_Emiss_hi1->Fill(Emiss,weight);
			h2p_Emiss_hi1_fine->Fill(Emiss,weight);
		      }
		    else
		      {
			h2p_Emiss_hi2->Fill(Emiss,weight);
			h2p_Emiss_hi2_fine->Fill(Emiss,weight);
		      }
		  }

		// Let's make a sanitized phi and sector
		double phi2_deg = vrec.Phi() * 180./M_PI;
		if (phi2_deg < -30.)
			phi2_deg += 360.;
		int sector2 = phi2_deg/60.;
		double theta2_deg = vrec.Theta() * 180./M_PI;

		h2p_phi2->Fill(phi2_deg,weight);
		h2p_theta2->Fill(theta2_deg,weight);
		h2p_theta2_bySec[sector2]->Fill(theta2_deg,weight);
		h2p_mom2->Fill(Pp_size[1],weight);

		// Histogram for the "apparent E*"
		double apparent_Estar = sqrt(sq(sqrt(sq(m_10B) + vcm.Mag2()) + Erec) -vlead.Mag2()) - m_11B;
		h2p_pmiss_appEstar->Fill(Pmiss_size[0],apparent_Estar,weight);
	}
	f1p->Close();
	f2p->Close();

	cerr << "The ep and epp integrals are: " << h1p_Pm->Integral() << " "  << h2p_Pm->Integral() << "\n";
	cerr << "Broken down by pmiss range...\n\n";
	for (int j=4 ; j<=28 ; j+=4)
	  {
	    double min=0.3 + 0.1*(j-4)/4.;
	    double max=0.3 + 0.1*j/4.;
	    cerr << min << " < pmiss < " << max << " : " << h1p_Pm->Integral(j-3,j) << " " << h2p_Pm->Integral(j-3,j) << "\n";
	  }

	// pp-to-p
	pp_to_p->BayesDivide(h2p_Pm,h1p_Pm);
	pp_to_p_coarse->BayesDivide(h2p_Pm_coarse,h1p_Pm_coarse);
	for (int binX=1 ; binX<=pp_to_p_2d->GetNbinsX() ; binX++)
	  for (int binY=1 ; binY<=pp_to_p_2d->GetNbinsY() ; binY++)
	    {
	      
	      double N_pp = h2p_pmiss_E1->GetBinContent(binX,binY);
	      double N_p =  h1p_pmiss_E1->GetBinContent(binX,binY);
	      double ratio = 0.;
	      if (N_p >0)
		ratio = N_pp / N_p;
	      pp_to_p_2d->SetBinContent(binX,binY,ratio);
	    }

	// Write out
	fo -> cd();
	pp_to_p->Write();
	pp_to_p_coarse->Write();
	pp_to_p_2d->Write();

	const double data_ep = 9175.;
	const double data_epp = 437.;
	const double pnorm = data_ep/h1p_Pm->Integral();
	const double ppnorm = pnorm;

	// Including a factor if we watn to rescale data to match epp luminosity
	const double renorm = data_epp/h2p_Pm->Integral()/(data_ep/h1p_Pm->Integral());
	TVectorT<double> renorm_factor(1);
	renorm_factor[0] = renorm;
	renorm_factor.Write("factor");

	// scale all the histograms, and write them out
	for (int i=0 ; i<h1p_list.size() ; i++)
	  {
	    h1p_list[i]->Scale(pnorm);
	    h1p_list[i]->Write();
	  }
	for (int i=0 ; i<h2p_list.size() ; i++)
	  {
	    h2p_list[i]->Scale(ppnorm);
	    h2p_list[i]->Write();
	  }

	fo->Close();

	return 0;
}

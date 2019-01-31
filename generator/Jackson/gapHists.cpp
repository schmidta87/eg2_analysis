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

	// Create histograms
	TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",40,1.,5.);
	h1p_list.push_back(h1p_QSq);
	TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",26,1.2,2.5);
	h1p_list.push_back(h1p_xB );
	TH1D * h1p_thetam = new TH1D("ep_thetam","ep;Theta_m [deg];Counts",40,10.,170.);
	h1p_list.push_back(h1p_thetam);
	TH1D * h1p_cthetam = new TH1D("ep_cthetam","ep;cos Theta_m;Counts",40,-1,1);
	h1p_list.push_back(h1p_cthetam);
	TH1D * h1p_Pmz =  new TH1D("ep_Pmz" ,"ep;pMiss_z [GeV];Counts",30,0.,0.4);
	h1p_list.push_back(h1p_Pmz);
	TH1D * h1p_PmT =  new TH1D("ep_PmT" ,"ep;pMiss_T [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmT);
	TH1D * h1p_Pmzq =  new TH1D("ep_Pmzq" ,"ep;pMiss_zq [GeV];Counts",30,-1.0,0.0);
	h1p_list.push_back(h1p_Pmzq);
	TH1D * h1p_PmTq =  new TH1D("ep_PmTq" ,"ep;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmTq);
	TH1D * h1p_Pm =  new TH1D("ep_Pm" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm);
	TH1D * h1p_Pm_10 =  new TH1D("ep_Pm_10" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_10);
	TH1D * h1p_Pm_15 =  new TH1D("ep_Pm_15" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_15);
	TH1D * h1p_Pm_17 =  new TH1D("ep_Pm_17" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_17);
	TH1D * h1p_Pm_20 =  new TH1D("ep_Pm_20" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_20);
	TH1D * h1p_cPeq = new TH1D("ep_cPeq","ep;cos(Theta_Peq);Counts",35,0.1,0.45);
	h1p_list.push_back(h1p_cPeq);
	TH1D * h1p_cPmq = new TH1D("ep_cPmq","ep;cos(Theta_Pmq);Counts",40,-1.,0.);
	h1p_list.push_back(h1p_cPmq);
	TH1D * h1p_cP1q = new TH1D("ep_cP1q","ep;cos(Theta_P1q);Counts",30,0.85,1.);
	h1p_list.push_back(h1p_cP1q);
	TH1D * h1p_cPm1 = new TH1D("ep_cPm1","ep;cos(Theta_Pm1);Counts",40,-1.,0.);
	h1p_list.push_back(h1p_cPm1);
	TH1D * h1p_theta1 = new TH1D("ep_theta1","ep;Theta_1 [deg];Counts",40,10.,90.);
	h1p_list.push_back(h1p_theta1);
	TH1D * h1p_thetae = new TH1D("ep_thetae","ep;Theta_e [deg];Counts",60,10.,40.);
	h1p_list.push_back(h1p_thetae);
	TH1D * h1p_mome = new TH1D("ep_mome","ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h1p_list.push_back(h1p_mome);
	TH1D * h1p_momez = new TH1D("ep_momez","ep;Mom_ez [GeV/c];Counts",60,2.7,4.7);
	h1p_list.push_back(h1p_momez);
	TH1D * h1p_momeT = new TH1D("ep_momeT","ep;Mom_eT [GeV/c];Counts",60,1.0,2.0);
	h1p_list.push_back(h1p_momeT);
	TH1D * h1p_momezq = new TH1D("ep_momezq","ep;Mom_ezq [GeV/c];Counts",60,0.0,2.0);
	h1p_list.push_back(h1p_momezq);
	TH1D * h1p_momeTq = new TH1D("ep_momeTq","ep;Mom_eTq [GeV/c];Counts",60,2.5,5.0);
	h1p_list.push_back(h1p_momeTq);
	TH1D * h1p_mom1 = new TH1D("ep_mom1","ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h1p_list.push_back(h1p_mom1);
	TH1D * h1p_mom1z = new TH1D("ep_mom1z","ep;Mom_1z [GeV/c];Counts",60,0.0,2.0);
	h1p_list.push_back(h1p_mom1z);
	TH1D * h1p_mom1T = new TH1D("ep_mom1T","ep;Mom_1T [GeV/c];Counts",60,0.0,2.0);
	h1p_list.push_back(h1p_mom1T);
	TH1D * h1p_mom1zq = new TH1D("ep_mom1zq","ep;Mom_1zq [GeV/c];Counts",60,0.5,3.0);
	h1p_list.push_back(h1p_mom1zq);
	TH1D * h1p_mom1Tq = new TH1D("ep_mom1Tq","ep;Mom_1Tq [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h1p_mom1Tq);
	TH1D * h1p_phi1 = new TH1D("ep_phi1","ep;Phi_1 [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phi1);
	TH1D * h1p_phie = new TH1D("ep_phie","ep;Phi_e [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phie);
	TH1D * h1p_q = new TH1D("ep_q","ep;q [GeV/c];Counts",40,1.0,3.0);
	h1p_list.push_back(h1p_q);
	TH1D * h1p_qz = new TH1D("ep_qz","ep;q_z [GeV/c];Counts",50,0.0,2.0);
	h1p_list.push_back(h1p_qz);
	TH1D * h1p_qT = new TH1D("ep_qT","ep;q_T [GeV/c];Counts",40,1.0,2.0);
	h1p_list.push_back(h1p_qT);
	TH1D * h1p_thetaq = new TH1D("ep_thetaq","ep;Theta_q [deg];Counts",50,30.,80.);
	h1p_list.push_back(h1p_thetaq);	
	TH1D * h2p_mom2 = new TH1D("epp_mom2","epp;Recoil Mom [GeV/c];Counts",17,0.35,1.2);
	h1p_list.push_back(h2p_mom2);
	TH1D * h2p_mom2z = new TH1D("epp_mom2z","ep;Mom_1z [GeV/c];Counts",60,0-1.2,1.2);
	h1p_list.push_back(h2p_mom2z);
	TH1D * h2p_mom2T = new TH1D("epp_mom2T","ep;Mom_1T [GeV/c];Counts",60,00,1.2);
	h1p_list.push_back(h2p_mom2T);
	TH1D * h2p_mom2zq = new TH1D("epp_mom2zq","ep;Mom_1zq [GeV/c];Counts",60,-1.2,1.2);
	h1p_list.push_back(h2p_mom2zq);
	TH1D * h2p_mom2Tq = new TH1D("epp_mom2Tq","ep;Mom_1Tq [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_mom2Tq);
	TH1D * h2p_pCM = new TH1D("epp_pCM","ep;p_CM [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_pCM);
	TH1D * h2p_prel = new TH1D("epp_prel","ep;p_rel [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_prel);
	
	TH2D * thetaq_vs_phiq = new TH2D("thetaq_vs_phiq","ep; theta_q [degrees]; phi_q [degrees]",40,30.,70.,60,-30.,330.);
	h1p_list.push_back(thetaq_vs_phiq);

	TH1D * h1p_thetae_bySec[6];
	TH1D * h1p_theta1_bySec[6];
	TH1D * h1p_thetaq_bySec[6];
	TH1D * h1p_mome_bySec[6];
	TH1D * h1p_q_bySec[6];
	TH1D * h1p_mom1_bySec[6];
	TH1D * h1p_cPmq_bySec[6];
	TH1D * h1p_Pm_bySec[6];
	TH1D * h1p_Pm_bySece[6];
	TH1D * h1p_Em_bySece[6];
	TH1D * h1p_Pm_bySece_tot[6];
	TH1D * h1p_Em_bySece_tot[6];
	TH1D * h1p_Pm_bySece_tot_10[6];
	TH1D * h1p_Pm_bySece_tot_15[6];
	TH1D * h1p_Pm_bySece_tot_17[6];
	TH1D * h1p_Pm_bySece_tot_20[6];

	for (int i=0 ; i<6 ; i++)
	{
		char temp[100];

		sprintf(temp,"ep_theta1_%d",i);
		h1p_theta1_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h1p_theta1_bySec[i]);

		sprintf(temp,"ep_thetae_%d",i);
		h1p_thetae_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,40.);
		h1p_list.push_back(h1p_thetae_bySec[i]);

		sprintf(temp,"ep_thetaq_%d",i);
		h1p_thetaq_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",40,30.,70.);
		h1p_list.push_back(h1p_thetaq_bySec[i]);

		sprintf(temp,"ep_mome_%d",i);
		h1p_mome_bySec[i] = new TH1D(temp,"ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
		h1p_list.push_back(h1p_mome_bySec[i]);

		sprintf(temp,"ep_mom1_%d",i);
		h1p_mom1_bySec[i] = new TH1D(temp,"ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
		h1p_list.push_back(h1p_mom1_bySec[i]);

		sprintf(temp,"ep_q_%d",i);
		h1p_q_bySec[i] = new TH1D(temp,"ep;q [GeV/c];Counts",40,1.0,3.0);
		h1p_list.push_back(h1p_q_bySec[i]);

		sprintf(temp,"ep_cPmq_%d",i);
		h1p_cPmq_bySec[i] = new TH1D(temp,"ep;cos(Theta_Pmq);Counts",40,-1.,0.);
		h1p_list.push_back(h1p_cPmq_bySec[i]);

		sprintf(temp,"ep_Pm_%d",i);
		h1p_Pm_bySec[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",10,0.35,0.55);
		h1p_list.push_back(h1p_Pm_bySec[i]);

		sprintf(temp,"ep_Pm_e_%d",i);
		h1p_Pm_bySece[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",10,0.35,0.55);
		h1p_list.push_back(h1p_Pm_bySece[i]);

		sprintf(temp,"ep_Em_e_%d",i);
		h1p_Em_bySece[i] = new TH1D(temp,"ep;EMiss [GeV];Counts",40,-0.2,0.6);
		h1p_list.push_back(h1p_Em_bySece[i]);

		sprintf(temp,"ep_Pm_etot_%d",i);
		h1p_Pm_bySece_tot[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot[i]);

		sprintf(temp,"ep_Em_etot_%d",i);
		h1p_Em_bySece_tot[i] = new TH1D(temp,"ep;EMiss [GeV];Counts",40,-0.2,0.6);
		h1p_list.push_back(h1p_Em_bySece_tot[i]);

		sprintf(temp,"ep_Pm_etot_10_%d",i);
		h1p_Pm_bySece_tot_10[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_10[i]);

		sprintf(temp,"ep_Pm_etot_15_%d",i);
		h1p_Pm_bySece_tot_15[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_15[i]);

		sprintf(temp,"ep_Pm_etot_17_%d",i);
		h1p_Pm_bySece_tot_17[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_17[i]);

		sprintf(temp,"ep_Pm_etot_20_%d",i);
		h1p_Pm_bySece_tot_20[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_20[i]);

	}

	// Now that all histograms have been defined, set them to Sumw2
	for (int i=0 ; i<h1p_list.size() ; i++)
	  h1p_list[i]->Sumw2();

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
	t1p->SetBranchAddress("q",q);
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

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		double omega = Q2/(2.*mN*Xb);

		// Let's make a sanitized phi and sector
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = phie_deg/60.;
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);

		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = phi1_deg/60.;
		double theta1_deg = vp.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 qu = vq.Unit();

		double phiq_deg = vq.Phi() * 180./M_PI;
		if (phiq_deg < -30.)
		  phiq_deg += 360.;
		double thetaq_deg = vq.Theta() * 180./M_PI;

		h1p_thetaq->Fill(thetaq_deg,weight);
		h1p_q->Fill(vq.Mag(),weight);
		h1p_qz ->Fill(vq.z(),weight);
		h1p_qT ->Fill(vq.XYvector().Mod(),weight);

		TVector3 vm = vp - vq;
		double thetam_deg = vm.Theta() * 180./M_PI;

		h1p_thetam->Fill(thetam_deg,weight);
		h1p_cthetam->Fill(cos(thetam_deg * M_PI/180.),weight);

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);

		h1p_mome->Fill(ve.Mag(),weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		h1p_momez ->Fill(ve.z(),weight);
		h1p_momeT ->Fill(ve.XYvector().Mod(),weight);
		h1p_mom1z ->Fill(vp.z(),weight);
		h1p_mom1T ->Fill(vp.XYvector().Mod(),weight);
		h1p_Pmz ->Fill(vm.z(),weight);
		h1p_PmT ->Fill(vm.XYvector().Mod(),weight);

		h1p_momezq ->Fill(ve.Dot(qu),weight);
		h1p_momeTq ->Fill(ve.Perp(qu),weight);
		h1p_mom1zq ->Fill(vp.Dot(qu),weight);
		h1p_mom1Tq ->Fill(vp.Perp(qu),weight);
		h1p_Pmzq ->Fill(vm.Dot(qu),weight);
		h1p_PmTq ->Fill(vm.Perp(qu),weight);

		h1p_cPeq->Fill(cos(ve.Angle(vq)),weight);
		h1p_cP1q->Fill(cos(vp.Angle(vq)),weight);
		h1p_cPm1->Fill(cos(vm.Angle(vp)),weight);

		if ((0.35<Pmiss_size[0]) and (Pmiss_size[0]<0.55))
		  {

		    h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		    h1p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		    h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		    h1p_mom1_bySec[sector]->Fill(Pp_size[0],weight);
		    h1p_thetaq_bySec[sec_e]->Fill(thetaq_deg,weight);
		    h1p_q_bySec[sec_e]->Fill(vq.Mag(),weight);
		    h1p_cPmq_bySec[sec_e]->Fill(cos(Pmiss_q_angle[0]*M_PI/180),weight); 
		    h1p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		    h1p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 
		    h1p_Em_bySece[sec_e]->Fill(Emiss,weight); 
		  }
		h1p_Pm_bySece_tot[sec_e]->Fill(Pmiss_size[0],weight); 
		h1p_Em_bySece_tot[sec_e]->Fill(Emiss,weight); 

		thetaq_vs_phiq->Fill(thetaq_deg,phiq_deg,weight);


		if (thetae_deg > 10)
		  {
		    h1p_Pm_10 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_10[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 15)
		  {
		    h1p_Pm_15 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_15[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 17.5)
		  {
		    h1p_Pm_17 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_17[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 20)
		  {
		    h1p_Pm_20 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_20[sec_e] ->Fill(Pmiss_size[0],weight);
		  }

	}

	f1p->Close();

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
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		double omega = Q2/(2.*mN*Xb);

		// Let's make a sanitized phi and sector
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = phie_deg/60.;
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);

		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = phi1_deg/60.;
		double theta1_deg = vp.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 qu = vq.Unit();

		double phiq_deg = vq.Phi() * 180./M_PI;
		if (phiq_deg < -30.)
		  phiq_deg += 360.;
		double thetaq_deg = vq.Theta() * 180./M_PI;

		h1p_thetaq->Fill(thetaq_deg,weight);
		h1p_q->Fill(vq.Mag(),weight);
		h1p_qz ->Fill(vq.z(),weight);
		h1p_qT ->Fill(vq.XYvector().Mod(),weight);

		TVector3 vm = vp - vq;
		double thetam_deg = vm.Theta() * 180./M_PI;

		h1p_thetam->Fill(thetam_deg,weight);
		h1p_cthetam->Fill(cos(thetam_deg * M_PI/180.),weight);

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);

		h1p_mome->Fill(ve.Mag(),weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		h1p_momez ->Fill(ve.z(),weight);
		h1p_momeT ->Fill(ve.XYvector().Mod(),weight);
		h1p_mom1z ->Fill(vp.z(),weight);
		h1p_mom1T ->Fill(vp.XYvector().Mod(),weight);
		h1p_Pmz ->Fill(vm.z(),weight);
		h1p_PmT ->Fill(vm.XYvector().Mod(),weight);

		h1p_momezq ->Fill(ve.Dot(qu),weight);
		h1p_momeTq ->Fill(ve.Perp(qu),weight);
		h1p_mom1zq ->Fill(vp.Dot(qu),weight);
		h1p_mom1Tq ->Fill(vp.Perp(qu),weight);
		h1p_Pmzq ->Fill(vm.Dot(qu),weight);
		h1p_PmTq ->Fill(vm.Perp(qu),weight);

		h1p_cPeq->Fill(cos(ve.Angle(vq)),weight);
		h1p_cP1q->Fill(cos(vp.Angle(vq)),weight);
		h1p_cPm1->Fill(cos(vm.Angle(vp)),weight);

		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
		h2p_mom2z ->Fill(vrec.z(),weight);
		h2p_mom2T ->Fill(vrec.XYvector().Mod(),weight);
		h2p_mom2zq ->Fill(vrec.Dot(qu),weight);
		h2p_mom2Tq ->Fill(vrec.Perp(qu),weight);
		
		TVector3 vcm = vm + vrec;
		TVector3 vrel = (vm - vrec)*(0.5);
		h2p_pCM ->Fill(vcm.Mag(),weight);
		h2p_prel ->Fill(vrel.Mag(),weight);
		

		if ((0.35<Pmiss_size[0]) and (Pmiss_size[0]<0.55))
		  {

		    h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		    h1p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		    h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		    h1p_mom1_bySec[sector]->Fill(Pp_size[0],weight);
		    h1p_thetaq_bySec[sec_e]->Fill(thetaq_deg,weight);
		    h1p_q_bySec[sec_e]->Fill(vq.Mag(),weight);
		    h1p_cPmq_bySec[sec_e]->Fill(cos(Pmiss_q_angle[0]*M_PI/180),weight); 
		    h1p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		    h1p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 
		    h1p_Em_bySece[sec_e]->Fill(Emiss,weight); 
		  }
		h1p_Pm_bySece_tot[sec_e]->Fill(Pmiss_size[0],weight); 
		h1p_Em_bySece_tot[sec_e]->Fill(Emiss,weight); 

		thetaq_vs_phiq->Fill(thetaq_deg,phiq_deg,weight);


		if (thetae_deg > 10)
		  {
		    h1p_Pm_10 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_10[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 15)
		  {
		    h1p_Pm_15 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_15[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 17.5)
		  {
		    h1p_Pm_17 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_17[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 20)
		  {
		    h1p_Pm_20 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_20[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		h2p_mom2->Fill(Pp_size[1],weight);

	}

	f2p->Close();

	cerr << "The ep integrals is: " << h1p_Pm->Integral() << "\n";
	cerr << "Broken down by pmiss range...\n\n";
	for (int j=4 ; j<=28 ; j+=4)
	  {
	    double min=0.3 + 0.1*(j-4)/4.;
	    double max=0.3 + 0.1*j/4.;
	    cerr << min << " < pmiss < " << max << " : " << h1p_Pm->Integral(j-3,j) << "\n";
	  }

	// Write out
	fo -> cd();

	const double data_ep = 9175.;
	const double data_epp = 437.;
	const double pnorm = data_ep/h1p_Pm->Integral();

	// scale all the histograms, and write them out
	for (int i=0 ; i<h1p_list.size() ; i++)
	  {
	    h1p_list[i]->Scale(pnorm);
	    h1p_list[i]->Write();
	  }

	fo->Close();

	return 0;
}

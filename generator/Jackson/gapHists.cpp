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
#include "AccMap.h"

using namespace std;

const double pmiss_lo=0.5;
const double pmiss_md=0.6;
const double pmiss_hi=0.7;

const double pmiss_cut=0.4;

const double acc_thresh=0.8;

int main(int argc, char ** argv)
{
	if (argc != 5)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   make_hists /path/to/1p/file /path/to/2p/file /path/to/map/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * f1p = new TFile(argv[1]);
	TFile * f2p = new TFile(argv[2]);
	AccMap proton_map(argv[3], "p");

	TFile * fo = new TFile(argv[4],"RECREATE");

	// Let's create a vector of all the histogram pointers so we can loop over them, save hassles
	vector<TH1*> h1p_list;

	// Create histograms
	TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",40,1.,5.);
	h1p_list.push_back(h1p_QSq);
	TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",26,1.2,2.5);
	h1p_list.push_back(h1p_xB );
	TH1D * h1p_thetam = new TH1D("ep_thetam","ep;Theta_m [deg];Counts",40,10.,170.);
	h1p_list.push_back(h1p_thetam);
	TH1D * h1p_Pm =  new TH1D("ep_Pm" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm);
	TH1D * h2p_Pm =  new TH1D("epp_Pm" ,"epp;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h2p_Pm);
	TH1D * h1p_theta1 = new TH1D("ep_theta1","ep;Theta_1 [deg];Counts",40,10.,90.);
	h1p_list.push_back(h1p_theta1);
	TH1D * h2p_theta1 = new TH1D("epp_theta1","epp;Theta_1 [deg];Counts",40,10.,90.);
	h1p_list.push_back(h2p_theta1);
	TH1D * h2p_theta2 = new TH1D("epp_theta2","epp;Theta_2 [deg];Counts",40,10.,90.);
	h1p_list.push_back(h2p_theta2);
	TH1D * h1p_thetae = new TH1D("ep_thetae","ep;Theta_e [deg];Counts",60,10.,40.);
	h1p_list.push_back(h1p_thetae);
	TH1D * h2p_thetae = new TH1D("epp_thetae","epp;Theta_e [deg];Counts",60,10.,40.);
	h1p_list.push_back(h2p_thetae);
	TH1D * h1p_mome = new TH1D("ep_mome","ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h1p_list.push_back(h1p_mome);
	TH1D * h1p_mom1 = new TH1D("ep_mom1","ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h1p_list.push_back(h1p_mom1);
	TH1D * h2p_mome = new TH1D("epp_mome","epp;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h1p_list.push_back(h2p_mome);
	TH1D * h2p_mom1 = new TH1D("epp_mom1","epp;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h1p_list.push_back(h2p_mom1);
	TH1D * h1p_phi1 = new TH1D("ep_phi1","ep;Phi_1 [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phi1);
	TH1D * h1p_phie = new TH1D("ep_phie","ep;Phi_e [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phie);
	TH1D * h1p_q = new TH1D("ep_q","ep;q [GeV/c];Counts",40,1.0,3.0);
	h1p_list.push_back(h1p_q);
	TH1D * h1p_thetaq = new TH1D("ep_thetaq","ep;Theta_q [deg];Counts",50,30.,80.);
	h1p_list.push_back(h1p_thetaq);	
	TH1D * h2p_mom2 = new TH1D("epp_mom2","epp;Recoil Mom [GeV/c];Counts",17,0.35,1.2);
	h1p_list.push_back(h2p_mom2);
	TH1D * h2p_prel = new TH1D("epp_prel","ep;p_rel [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_prel);

	TH1D * h2p_e1 =  new TH1D("epp_e1" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1);
	TH1D * h2p_e1_lo =  new TH1D("epp_e1_lo" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1_lo);
	TH1D * h2p_e1_md =  new TH1D("epp_e1_md" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1_md);
	TH1D * h2p_e1_hi =  new TH1D("epp_e1_hi" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1_hi);

	
	TH2D * thetaq_vs_phiq = new TH2D("thetaq_vs_phiq","ep; theta_q [degrees]; phi_q [degrees]",40,30.,70.,60,-30.,330.);
	h1p_list.push_back(thetaq_vs_phiq);
	TH2D * theta1_vs_phi1 = new TH2D("theta1_vs_phi1","ep; theta_1 [degrees]; phi_1 [degrees]",80,10.,90.,60,-30.,330.);
	h1p_list.push_back(theta1_vs_phi1);
	TH2D * thetae_vs_phie = new TH2D("thetae_vs_phie","ep; theta_e [degrees]; phi_e [degrees]",60,10.,40.,60,-30.,330.);
	h1p_list.push_back(thetae_vs_phie);

	TH2D * prel_vs_Qsq = new TH2D("prel_vs_Qsq","ep; p_rel [GeV/c]; QSq [GeV^2]",40,0.0,1.0,40,1.,5.);
	h1p_list.push_back(prel_vs_Qsq);
	TH2D * prel_vs_xB = new TH2D("prel_vs_xB","ep; p_rel [GeV/c]; xB",40,0.0,1.0,26,1.2,2.5);
	h1p_list.push_back(prel_vs_xB);
	TH2D * prel_vs_Pm = new TH2D("prel_vs_Pm","ep; p_rel [GeV/c]; Pm [GeV/c]",40,0.0,1.0,28,0.3,1.0);
	h1p_list.push_back(prel_vs_Pm);
	TH2D * prel_vs_Em = new TH2D("prel_vs_Em","ep; p_rel [GeV/c]; Em [GeV]",40,0.0,1.0,40,-0.2,0.6);
	h1p_list.push_back(prel_vs_Em);
	TH2D * prel_vs_mom2 = new TH2D("prel_vs_mom2","ep; p_rel [GeV/c]; mom2 [GeV/c]",40,0.0,1.0,17,0.35,1.2);
	h1p_list.push_back(prel_vs_mom2);

	TH2D * prec_vs_e1 = new TH2D("prec_vs_e1","ep; p_rec [GeV/c]; e1 [GeV]",17,0.35,1.2,40,0.6,1.6);
	h1p_list.push_back(prec_vs_e1);

	TH2D * Pm_vs_Thetam = new TH2D("Pm_vs_Thetam","ep; p_miss [GeV/c]; Theta_m [degrees]",28,0.3,1.0,40,10.,170.);
	h1p_list.push_back(Pm_vs_Thetam);
	TH2D * Pmz_vs_PmT = new TH2D("Pmz_vs_PmT","ep; p_mz [GeV/c]; p_mT [GeV/c]",30,0.,0.4,30,0.0,1.0);
	h1p_list.push_back(Pmz_vs_PmT);
	TH2D * Pmzq_vs_PmTq = new TH2D("Pmzq_vs_PmTq","ep; p_mzq [GeV/c]; p_mTq [GeV/c]",30,-1.0,0.0,30,0.0,1.0);
	h1p_list.push_back(Pmzq_vs_PmTq);

	TH2D * P1_vs_Theta1 = new TH2D("P1_vs_Theta1","ep; p_1 [GeV/c]; Theta_1 [degrees]",40,0.4,2.4,40,10.,90.);
	h1p_list.push_back(P1_vs_Theta1);

	TH1D * h1p_thetae_bySec[6];
	TH1D * h1p_theta1_bySec[6];
	TH1D * h1p_mome_bySec[6];
	TH1D * h1p_mom1_bySec[6];
	TH1D * h1p_Pm_bySec[6];
	TH1D * h1p_Pm_bySece[6];
	TH1D * h2p_thetae_bySec[6];
	TH1D * h2p_theta1_bySec[6];
	TH1D * h2p_theta2_bySec[6];
	TH1D * h2p_mome_bySec[6];
	TH1D * h2p_mom1_bySec[6];
	TH1D * h2p_mom2_bySec[6];
	TH1D * h2p_Pm_bySec[6];
	TH1D * h2p_Pm_bySece[6];
	TH1D * h1p_Em_bySece[6];
	TH1D * h1p_q_bySece[6];
	TH1D * h1p_thetaq_bySec[6];

	TH2D * P1_vs_Theta1_bySec[6];


	for (int i=0 ; i<6 ; i++)
	{
		char temp[100];

		sprintf(temp,"ep_theta1_%d",i);
		h1p_theta1_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h1p_theta1_bySec[i]);

		sprintf(temp,"ep_thetae_%d",i);
		h1p_thetae_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,40.);
		h1p_list.push_back(h1p_thetae_bySec[i]);

		sprintf(temp,"ep_mome_%d",i);
		h1p_mome_bySec[i] = new TH1D(temp,"ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
		h1p_list.push_back(h1p_mome_bySec[i]);

		sprintf(temp,"ep_mom1_%d",i);
		h1p_mom1_bySec[i] = new TH1D(temp,"ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
		h1p_list.push_back(h1p_mom1_bySec[i]);

		sprintf(temp,"ep_Pm_%d",i);
		h1p_Pm_bySec[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySec[i]);

		sprintf(temp,"ep_Pm_e_%d",i);
		h1p_Pm_bySece[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece[i]);

		sprintf(temp,"epp_theta1_%d",i);
		h2p_theta1_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h2p_theta1_bySec[i]);

		sprintf(temp,"epp_theta2_%d",i);
		h2p_theta2_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h2p_theta2_bySec[i]);

		sprintf(temp,"epp_thetae_%d",i);
		h2p_thetae_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,40.);
		h1p_list.push_back(h2p_thetae_bySec[i]);

		sprintf(temp,"epp_mome_%d",i);
		h2p_mome_bySec[i] = new TH1D(temp,"epp;Mom_e [GeV/c];Counts",40,3.0,5.0);
		h1p_list.push_back(h2p_mome_bySec[i]);

		sprintf(temp,"epp_mom1_%d",i);
		h2p_mom1_bySec[i] = new TH1D(temp,"epp;Mom_1 [GeV/c];Counts",40,0.4,2.4);
		h1p_list.push_back(h2p_mom1_bySec[i]);

		sprintf(temp,"epp_mom2_%d",i);
		h2p_mom2_bySec[i] = new TH1D(temp,"epp;Mom_2 [GeV/c];Counts",40,0.35,1.2);
		h1p_list.push_back(h2p_mom2_bySec[i]);

		sprintf(temp,"epp_Pm_%d",i);
		h2p_Pm_bySec[i] = new TH1D(temp,"epp;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h2p_Pm_bySec[i]);

		sprintf(temp,"epp_Pm_e_%d",i);
		h2p_Pm_bySece[i] = new TH1D(temp,"epp;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h2p_Pm_bySece[i]);

		sprintf(temp,"ep_Em_e_%d",i);
		h1p_Em_bySece[i] = new TH1D(temp,"ep;EMiss [GeV];Counts",40,-0.2,0.6);
		h1p_list.push_back(h1p_Em_bySece[i]);

		sprintf(temp,"ep_q_%d",i);
		h1p_q_bySece[i] = new TH1D(temp,"ep;q [GeV/c];Counts",40,1.0,3.0);
		h1p_list.push_back(h1p_q_bySece[i]);

		sprintf(temp,"ep_thetaq_%d",i);
		h1p_thetaq_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",40,30.,70.);
		h1p_list.push_back(h1p_thetaq_bySec[i]);

		sprintf(temp,"p1_vs_theta1_%d",i);
		P1_vs_Theta1_bySec[i] = new TH2D(temp,"ep; p_1 [GeV/c]; Theta_1 [degrees]",40,0.4,2.4,40,10.,90.);
		h1p_list.push_back(P1_vs_Theta1_bySec[i]);

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
		if (Pmiss_size[0] < pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;

		// Apply an additional fiducial cut that the map acc must be > acc_thresh
		if ( proton_map.accept(vp) < acc_thresh)
		  continue;

		// Let's make a sanitized phi and sector
		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vp.Theta() * 180./M_PI;

		if ((sector==90) and (theta1_deg < 52))
		  continue;
		if ((sector==92) and ((theta1_deg < 45) or (theta1_deg > 72)))
		  continue;
		if ((sector==93) and (theta1_deg < 50))
		  continue;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);

		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = clas_sector(phie_deg);
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 qu = vq.Unit();

		double phiq_deg = vq.Phi() * 180./M_PI;
		if (phiq_deg < -30.)
		  phiq_deg += 360.;
		double thetaq_deg = vq.Theta() * 180./M_PI;

		h1p_thetaq->Fill(thetaq_deg,weight);
		h1p_q->Fill(vq.Mag(),weight);

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);

		double omega = Q2/(2.*mN*Xb);

		TVector3 vm = vp - vq;
		double thetam_deg = vm.Theta() * 180./M_PI;

		h1p_thetam->Fill(thetam_deg,weight);

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);

		h1p_mome->Fill(ve.Mag(),weight);
		h1p_mom1->Fill(Pp_size[0],weight);

		Pm_vs_Thetam->Fill(Pmiss_size[0],thetam_deg,weight);
		
		P1_vs_Theta1->Fill(Pp_size[0],theta1_deg,weight);
		P1_vs_Theta1_bySec[sector]->Fill(Pp_size[0],theta1_deg,weight);
	
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1_bySec[sector]->Fill(Pp_size[0],weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		h1p_thetaq_bySec[sec_e]->Fill(thetaq_deg,weight);
		h1p_q_bySece[sec_e]->Fill(vq.Mag(),weight);
		h1p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		h1p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 
		h1p_Em_bySece[sec_e]->Fill(Emiss,weight); 	

		thetaq_vs_phiq->Fill(thetaq_deg,phiq_deg,weight);
		theta1_vs_phi1->Fill(theta1_deg,phi1_deg,weight);
		thetae_vs_phie->Fill(thetae_deg,phie_deg,weight);

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
		if (Pmiss_size[0] < pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;

		// Apply an additional fiducial cut that the map acc must be > acc_thresh
		if ( proton_map.accept(vp) < acc_thresh)
		  continue;

		// Let's make a sanitized phi and sector
		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vp.Theta() * 180./M_PI;

		if ((sector==90) and (theta1_deg < 52))
		  continue;
		if ((sector==92) and ((theta1_deg < 45) or (theta1_deg > 72)))
		  continue;
		if ((sector==93) and (theta1_deg < 50))
		  continue;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);

		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = clas_sector(phie_deg);
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 qu = vq.Unit();

		double phiq_deg = vq.Phi() * 180./M_PI;
		if (phiq_deg < -30.)
		  phiq_deg += 360.;
		double thetaq_deg = vq.Theta() * 180./M_PI;

		h1p_thetaq->Fill(thetaq_deg,weight);
		h1p_q->Fill(vq.Mag(),weight);

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);

		double omega = Q2/(2.*mN*Xb);
		double e1 = sqrt(vp.Mag2()+mN*mN)-omega;

		TVector3 vm = vp - vq;
		TVector3 mu = vm.Unit();
		double thetam_deg = vm.Theta() * 180./M_PI;

		h1p_thetam->Fill(thetam_deg,weight);

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);

		h1p_mome->Fill(ve.Mag(),weight);
		h1p_mom1->Fill(Pp_size[0],weight);

		Pm_vs_Thetam->Fill(Pmiss_size[0],thetam_deg,weight);
		
		P1_vs_Theta1->Fill(Pp_size[0],theta1_deg,weight);
		P1_vs_Theta1_bySec[sector]->Fill(Pp_size[0],theta1_deg,weight);

		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1_bySec[sector]->Fill(Pp_size[0],weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		h1p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		h1p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 
		h1p_Em_bySece[sec_e]->Fill(Emiss,weight); 

		h1p_thetaq_bySec[sec_e]->Fill(thetaq_deg,weight);
		h1p_q_bySece[sec_e]->Fill(vq.Mag(),weight);

		thetaq_vs_phiq->Fill(thetaq_deg,phiq_deg,weight);
		theta1_vs_phi1->Fill(theta1_deg,phi1_deg,weight);
		thetae_vs_phie->Fill(thetae_deg,phie_deg,weight);

		if (Pp_size[1] < 0.35)
		  continue;

		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
		double phi2_deg = vrec.Phi() * 180./M_PI;
		if (phi2_deg < -30.)
		  phi2_deg += 360.;
		int sector2 = clas_sector(phi2_deg);
		double theta2_deg = vrec.Theta() * 180./M_PI;

		if(true){
		  if (!accept_proton(vrec))
		    continue;

		  // Apply an additional fiducial cut that the map acc must be > acc_thresh
		  if ( proton_map.accept(vrec) < acc_thresh)
		    continue;
		}

		h2p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h2p_mom1_bySec[sector]->Fill(Pp_size[0],weight);
		h2p_theta2_bySec[sector2]->Fill(theta2_deg,weight);
		h2p_mom2_bySec[sector2]->Fill(vrec.Mag(),weight);
		h2p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h2p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		h2p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		h2p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 


		h2p_theta1->Fill(theta1_deg,weight);
		h2p_mom1->Fill(Pp_size[0],weight);
		h2p_theta2->Fill(theta2_deg,weight);
		h2p_mom2->Fill(vrec.Mag(),weight);
		h2p_thetae->Fill(thetae_deg,weight);
		h2p_mome->Fill(ve.Mag(),weight);
		h2p_Pm->Fill(Pmiss_size[0],weight); 


		h2p_e1 ->Fill(e1,weight);
		if (Pp_size[1] < 0.45)
		  h2p_e1_lo ->Fill(e1,weight);
		else if (Pp_size[1] < 0.55)
		  h2p_e1_md ->Fill(e1,weight);
		else
		  h2p_e1_hi ->Fill(e1,weight);


		TVector3 vcm = vm + vrec;
		TVector3 ucm = vcm.Unit();
		TVector3 vrel = (vm - vrec)*(0.5);
		TVector3 urel = vrel.Unit();
		Double_t prel = vrel.Mag();
		
		prel_vs_Qsq->Fill(prel,Q2,weight);
		prel_vs_xB->Fill(prel,Xb,weight);
		prel_vs_Pm->Fill(prel,Pmiss_size[0],weight);
		prel_vs_Em->Fill(prel,Emiss,weight);
		prel_vs_mom2->Fill(prel,Pp_size[1],weight);

		prec_vs_e1->Fill(Pp_size[1],e1,weight);

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

	const double data_ep = 5561.;
	const double data_epp = 359.;
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

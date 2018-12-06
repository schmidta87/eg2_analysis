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

	// Create histograms
	TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",40,1.,5.);
	h1p_QSq->Sumw2();
	TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",26,1.2,2.5);
	h1p_xB ->Sumw2();
	TH1D * h1p_Pm =  new TH1D("ep_Pm" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_Pm ->Sumw2();
	TH1D * h1p_Pmq = new TH1D("ep_Pmq","ep;Theta_Pmq [deg];Counts",40,100.,180.);
	h1p_Pmq->Sumw2();
	TH1D * h1p_phi1 = new TH1D("ep_phi1","ep;Phi [deg];Counts",60,-30.,330.);
	h1p_phi1->Sumw2();
	TH1D * h1p_phie = new TH1D("ep_phie","ep;Phi [deg];Counts",60,-30.,330.);
	h1p_phie->Sumw2();
	TH1D * h1p_theta1 = new TH1D("ep_theta1","ep;Theta [deg];Counts",60,10.,130.);
	h1p_theta1->Sumw2();
	TH1D * h1p_thetae = new TH1D("ep_thetae","ep;Theta [deg];Counts",60,10.,40.);
	h1p_thetae->Sumw2();
	TH1D * h1p_Emiss = new TH1D("ep_Emiss","ep;Emiss [GeV];Counts",50,0.,2);
	h1p_Emiss->Sumw2();
	TH1D * h2p_QSq = new TH1D("epp_QSq","epp;QSq [GeV^2];Counts",40,1.,5.);
	h2p_QSq->Sumw2();
	TH1D * h2p_xB =  new TH1D("epp_xB" ,"epp;xB;Counts",26,1.2,2.5);
	h2p_xB ->Sumw2();
	TH1D * h2p_Pm =  new TH1D("epp_Pm" ,"epp;pMiss [GeV];Counts",28,0.3,1.0);
	h2p_Pm ->Sumw2();
	TH1D * h2p_Pmq = new TH1D("epp_Pmq","epp;Theta_Pmq [deg];Counts",40,100.,180.);
	h2p_Pmq->Sumw2();
	TH1D * h2p_Pmr = new TH1D("epp_Pmr","epp;Theta_Pmr [deg];Counts",40,100.,180.);
	h2p_Pmr->Sumw2();
	TH1D * h2p_Pr = new TH1D("epp_Pr","epp;pRec [GeV];Counts",26,0.35,1.0);
	h2p_Pr->Sumw2();
	TH1D * h2p_phi1 = new TH1D("epp_phi1","ep;Phi [deg];Counts",60,-30.,330.);
	h2p_phi1->Sumw2();
	TH1D * h2p_phi2 = new TH1D("epp_phi2","ep;Phi [deg];Counts",60,-30.,330.);
	h2p_phi2->Sumw2();
	TH1D * h2p_phie = new TH1D("epp_phie","ep;Phi [deg];Counts",60,-30.,330.);
	h2p_phie->Sumw2();
	TH1D * h2p_thetae = new TH1D("epp_thetae","ep;Theta [deg];Counts",60,10.,40.);
	h2p_thetae->Sumw2();
	TH1D * h2p_theta1 = new TH1D("epp_theta1","ep;Theta [deg];Counts",60,10.,130.);
	h2p_theta1->Sumw2();
	TH1D * h2p_theta2 = new TH1D("epp_theta2","ep;Theta [deg];Counts",60,10.,130.);
	h2p_theta2->Sumw2();
	TH1D * h2p_Emiss = new TH1D("epp_Emiss","epp;Emiss [GeV];Counts",50,0.,2);
	h2p_Emiss->Sumw2();

	TH1D * h1p_thetae_bySec[6];
	TH1D * h1p_theta1_bySec[6];
	TH1D * h2p_theta1_bySec[6];
	TH1D * h2p_theta2_bySec[6];
	for (int i=0 ; i<6 ; i++)
	{
		char temp[100];

		sprintf(temp,"ep_theta1_%d",i);
		h1p_theta1_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,130.);
		h1p_theta1_bySec[i]->Sumw2();

		sprintf(temp,"ep_thetae_%d",i);
		h1p_thetae_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,40.);
		h1p_thetae_bySec[i]->Sumw2();

		sprintf(temp,"epp_theta1_%d",i);
		h2p_theta1_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h2p_theta1_bySec[i]->Sumw2();

		sprintf(temp,"epp_theta2_%d",i);
		h2p_theta2_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h2p_theta2_bySec[i]->Sumw2();
	}

	// pp2p hist
	TGraphAsymmErrors * pp_to_p = new TGraphAsymmErrors();
	pp_to_p->SetName("pp_to_p");
	pp_to_p->SetTitle("pp_to_p;p_miss [GeV];pp_to_p ratio");

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

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);

		// Let's make a sanitized phi and sector
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = phie_deg/60.;
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);

		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = phi1_deg/60.;
		double theta1_deg = vp.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);

		// Let's figure out missing energy! 
		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pmiss_size[0]*Pmiss_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);
		h1p_Emiss->Fill(Emiss,weight);
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

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);

		// Let's make a sanitized phi and sector
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = phie_deg/60.;
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);


		double phi1_deg = vlead.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = phi1_deg/60.;
		double theta1_deg = vlead.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);

		// Let's figure out missing energy! 
		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pmiss_size[0]*Pmiss_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);
		h1p_Emiss->Fill(Emiss,weight);

		// Make a check on the recoils
		if (fabs(Rp[1][2]+22.25)>2.25)
			continue;
		if (Pp_size[1] < 0.35)
			continue;
		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);      
		if (!accept_proton(vrec))
			continue;

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 vmiss = vlead - vq;

		h2p_QSq->Fill(Q2,weight);
		h2p_xB ->Fill(Xb,weight);
		h2p_Pm ->Fill(Pmiss_size[0],weight);
		h2p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h2p_Pr ->Fill(Pp_size[1],weight);
		h2p_Pmr->Fill(vmiss.Angle(vrec)*180./M_PI,weight);

		h2p_phie->Fill(phie_deg,weight);
		h2p_thetae->Fill(thetae_deg,weight);

		h2p_phi1->Fill(phi1_deg,weight);
		h2p_theta1->Fill(theta1_deg,weight);
		h2p_theta1_bySec[sector]->Fill(theta1_deg,weight);

		h2p_Emiss->Fill(Emiss,weight);

		// Let's make a sanitized phi and sector
		double phi2_deg = vrec.Phi() * 180./M_PI;
		if (phi2_deg < -30.)
			phi2_deg += 360.;
		int sector2 = phi2_deg/60.;
		double theta2_deg = vrec.Theta() * 180./M_PI;

		h2p_phi2->Fill(phi2_deg,weight);
		h2p_theta2->Fill(theta2_deg,weight);
		h2p_theta2_bySec[sector2]->Fill(theta2_deg,weight);

	}
	f1p->Close();
	f2p->Close();

	// pp-to-p
	pp_to_p->BayesDivide(h2p_Pm,h1p_Pm);

	// Write out
	fo -> cd();
	pp_to_p->Write();

	double pnorm = 9170;
	double ppnorm = 437;
						
	h1p_QSq->    Scale(pnorm/h1p_QSq->     Integral());
	h1p_xB->     Scale(pnorm/h1p_xB->      Integral());
	h1p_Pm->     Scale(pnorm/h1p_Pm->      Integral());
	h1p_Pmq->    Scale(pnorm/h1p_Pmq->     Integral());
	h1p_phi1->   Scale(pnorm/h1p_phi1->    Integral());
	h1p_theta1-> Scale(pnorm/h1p_theta1->  Integral());
	h1p_phie->   Scale(pnorm/h1p_phie->    Integral());
	h1p_thetae-> Scale(pnorm/h1p_thetae->  Integral());
	h1p_Emiss->  Scale(pnorm/h1p_Emiss->   Integral());
	h2p_QSq->    Scale(ppnorm/h2p_QSq->    Integral());
	h2p_xB->     Scale(ppnorm/h2p_xB->     Integral());
	h2p_Pm->     Scale(ppnorm/h2p_Pm->     Integral());
	h2p_Pmq->    Scale(ppnorm/h2p_Pmq->    Integral());
	h2p_Pr->     Scale(ppnorm/h2p_Pr->     Integral());
	h2p_Pmr->    Scale(ppnorm/h2p_Pmr->    Integral());
	h2p_phi1->   Scale(ppnorm/h2p_phi1->   Integral());
	h2p_theta1-> Scale(ppnorm/h2p_theta1-> Integral());
	h2p_phi2->   Scale(ppnorm/h2p_phi2->   Integral());
	h2p_theta2-> Scale(ppnorm/h2p_theta2-> Integral());
	h2p_phie->   Scale(ppnorm/h2p_phie->   Integral());
	h2p_thetae-> Scale(ppnorm/h2p_thetae-> Integral());
	h2p_Emiss->  Scale(ppnorm/h2p_Emiss->  Integral());

	h1p_QSq->Write();
	h1p_xB ->Write();
	h1p_Pm ->Write();
	h1p_Pmq->Write();
	h1p_phi1->Write();
	h1p_theta1->Write();
	h1p_phie->Write();
	h1p_thetae->Write();
	h1p_Emiss->Write();
	h2p_QSq->Write();
	h2p_xB ->Write();
	h2p_Pm ->Write();
	h2p_Pmq->Write();
	h2p_Pr ->Write();
	h2p_Pmr->Write();
	h2p_phi1->Write();
	h2p_theta1->Write();
	h2p_phi2->Write();
	h2p_theta2->Write();
	h2p_phie->Write();
	h2p_thetae->Write();
	h2p_Emiss->Write();

	for (int i=0 ; i<6 ; i++)
	{
		
		h1p_theta1_bySec[i]->Scale(ppnorm/h1p_theta1_bySec[i]->Integral());
		h1p_thetae_bySec[i]->Scale(ppnorm/h1p_thetae_bySec[i]->Integral());
		h2p_theta1_bySec[i]->Scale(ppnorm/h2p_theta1_bySec[i]->Integral());	  
		h2p_theta2_bySec[i]->Scale(ppnorm/h2p_theta2_bySec[i]->Integral());
		h1p_theta1_bySec[i]->Write();
		h1p_thetae_bySec[i]->Write();
		h2p_theta1_bySec[i]->Write();
		h2p_theta2_bySec[i]->Write();
	}

	fo->Close();

	return 0;
}

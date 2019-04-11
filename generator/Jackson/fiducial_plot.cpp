#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TVectorT.h"

#include "Nuclear_Info.h"
#include "fiducials.h"
#include "helpers.h"
#include "AccMap.h"

using namespace std;

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
	TFile * fo = new TFile(argv[4],"RECREATE");

	// We'll need to get acceptance maps in order to do a fiducial cut on minimum acceptance
	AccMap proton_map(argv[3], "p");

	TGraph * lead_fail = new TGraph();
	lead_fail->SetNameTitle("lead_fail","lead_fail");
	TGraph * lead_pass = new TGraph();
	lead_pass->SetNameTitle("lead_pass","lead_pass");
	TGraph * rec_fail = new TGraph();
	rec_fail->SetNameTitle("rec_fail","rec_fail");
	TGraph * rec_pass = new TGraph();
	rec_pass->SetNameTitle("rec_pass","rec_pass");

	TGraph * ep_theta_p_fail = new TGraph();
	ep_theta_p_fail->SetNameTitle("ep_theta_p_fail","ep_theta_p_fail");
	TGraph * ep_theta_p_pass = new TGraph();
	ep_theta_p_pass->SetNameTitle("ep_theta_p_pass","ep_theta_p_pass");
	TGraph * epp_theta_p_fail = new TGraph();
	epp_theta_p_fail->SetNameTitle("epp_theta_p_fail","epp_theta_p_fail");
	TGraph * epp_theta_p_pass = new TGraph();
	epp_theta_p_pass->SetNameTitle("epp_theta_p_pass","epp_theta_p_pass");
	TGraph * epp_theta_p_fail_rec = new TGraph();
	epp_theta_p_fail_rec->SetNameTitle("epp_theta_p_fail_rec","epp_theta_p_fail_rec");
	TGraph * epp_theta_p_pass_rec = new TGraph();
	epp_theta_p_pass_rec->SetNameTitle("epp_theta_p_pass_rec","epp_theta_p_pass_rec");

	TH1D * ep_pmiss = new TH1D("ep_pmiss","ep;pmiss;Counts",6,0.4,1.0);
	TH1D * epp_pmiss = new TH1D("epp_pmiss","epp;pmiss;Counts",6,0.4,1.0);

	ep_pmiss->Sumw2();
	epp_pmiss->Sumw2();

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
		if (Pmiss_size[0]<pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		if (!accept_electron(ve))
		  continue;

		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		double theta1_deg = vp.Theta() * 180./M_PI;
		if (!accept_proton(vp))
		  {
		    lead_fail->SetPoint(lead_fail->GetN(),theta1_deg,phi1_deg);
		    ep_theta_p_fail->SetPoint(ep_theta_p_fail->GetN(),theta1_deg,vp.Mag());
		    continue;
		  }
		lead_pass->SetPoint(lead_pass->GetN(),theta1_deg,phi1_deg);
		ep_theta_p_pass->SetPoint(ep_theta_p_pass->GetN(),theta1_deg,vp.Mag());

		if ( proton_map.accept(vp) > acc_thresh )
		  ep_pmiss->Fill(Pmiss_size[0]);

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
		if (Pmiss_size[0]<pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		if (!accept_electron(ve))
		  continue;

		TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
		double phi1_deg = vlead.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		double theta1_deg = vlead.Theta() * 180./M_PI;
		if (!accept_proton(vlead))
		  {
		    lead_fail->SetPoint(lead_fail->GetN(),theta1_deg,phi1_deg);
		    epp_theta_p_fail->SetPoint(epp_theta_p_fail->GetN(),theta1_deg,vlead.Mag());
		    continue;
		  }
		lead_pass->SetPoint(lead_pass->GetN(),theta1_deg,phi1_deg);
		epp_theta_p_pass->SetPoint(epp_theta_p_pass->GetN(),theta1_deg,vlead.Mag());

		if ( proton_map.accept(vlead) > acc_thresh )
		  {
		    ep_pmiss->Fill(Pmiss_size[0]);
		  }

		// Make a check on the recoils
		if (fabs(Rp[1][2]+22.25)>2.25)
		  continue;
		if (Pp_size[1] < 0.35)
		  continue;

		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);      
		double phi2_deg = vrec.Phi() * 180./M_PI;
		if (phi2_deg < -30.)
			phi2_deg += 360.;
		int sector2 = clas_sector(phi2_deg);
		double theta2_deg = vrec.Theta() * 180./M_PI;

		if (!accept_proton(vrec))
		  {
		    rec_fail->SetPoint(rec_fail->GetN(),theta2_deg,phi2_deg);
		    epp_theta_p_fail_rec->SetPoint(epp_theta_p_fail_rec->GetN(),theta2_deg,vrec.Mag());
		    continue;
		  }
		rec_pass->SetPoint(rec_pass->GetN(),theta2_deg,phi2_deg);
		epp_theta_p_pass_rec->SetPoint(epp_theta_p_pass_rec->GetN(),theta2_deg,vrec.Mag());
		
		if ( proton_map.accept(vlead) > acc_thresh and proton_map.accept(vrec) > acc_thresh )
		  {
		    epp_pmiss->Fill(Pmiss_size[0]);
		  }

	}
	f1p->Close();
	f2p->Close();


	// Write out
	fo -> cd();
	
	lead_fail->Write();
	lead_pass->Write();
	rec_fail->Write();
	rec_pass->Write();
	ep_theta_p_fail->Write();
	ep_theta_p_pass->Write();
	epp_theta_p_fail->Write();
	epp_theta_p_pass->Write();
	epp_theta_p_fail_rec->Write();
	epp_theta_p_pass_rec->Write();

	ep_pmiss->Write();
	epp_pmiss->Write();

	fo->Close();
	cerr << "The ep and epp integrals are: " << ep_pmiss->Integral() << " "  << epp_pmiss->Integral() << "\n";
	
	return 0;
}

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

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

    }
  f1p->Close();
  f2p->Close();

  // Write out
  fo -> cd();
  h1p_QSq->Write();
  h1p_xB ->Write();
  h1p_Pm ->Write();
  h1p_Pmq->Write();
  h2p_QSq->Write();
  h2p_xB ->Write();
  h2p_Pm ->Write();
  h2p_Pmq->Write();
  h2p_Pr ->Write();
  h2p_Pmr->Write();
  fo->Close();
  
  return 0;
}

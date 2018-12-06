#include <cstdlib>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

#include "AccMap.h"
#include "fiducials.h"
#include "Nuclear_Info.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 5)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tsimulator /path/to/gen/file /path/to/map/file /path/to/output/file [# of protons (1 or 2)]\n\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  AccMap pMap(argv[2]);
  AccMap eMap(argv[2],"e");
  TFile * outfile = new TFile(argv[3],"RECREATE");

  // Histograms
  TH1D * h_rec_p_all = new TH1D("rec_p_all","All recoil protons;pMiss [GeV];Counts",28,0.3,1.0);
  h_rec_p_all->Sumw2();
  TH1D * h_rec_p_acc = new TH1D("rec_p_acc","Accepted recoil protons;pMiss [GeV];Counts",28,0.3,1.0);
  h_rec_p_acc->Sumw2();
  TGraphAsymmErrors * rec_p_rat = new TGraphAsymmErrors();
  rec_p_rat->SetName("rec_p_rat");
  rec_p_rat->SetTitle("Acceptance Ratio;p_miss [GeV];Acceptance Ratio");

  Int_t nmb=atoi(argv[4]);
  if ((nmb < 1) || (nmb>2))
    {
      cerr << "The number of protons specified must be 1 or 2. You entered: " << nmb << "\n"
	   << "Try again.\n\n";
      return -2;
    }

  // Input Tree
  TTree * inTree = (TTree*)infile->Get("genT");
  Double_t gen_pe[3], gen_pLead[3], gen_pRec[3], gen_weight;
  Int_t lead_type, rec_type;
  inTree->SetBranchAddress("lead_type",&lead_type);
  inTree->SetBranchAddress("rec_type",&rec_type);
  inTree->SetBranchAddress("weight",&gen_weight);
  inTree->SetBranchAddress("pe",gen_pe);
  inTree->SetBranchAddress("pLead",gen_pLead);
  inTree->SetBranchAddress("pRec",gen_pRec);

  // Output Tree
  TTree * T = new TTree("T","Simulated Data Tree");
  Float_t Q2, Xb, Pe[3], Pe_size, theta_e, phi_e, Pp[2][3], Pp_size[2], pq_angle[2], Ep[2], theta_p[2], phi_p[2], nu, q[3];
  Float_t Pmiss_q_angle[2], Pmiss_size[2], Pmiss[2][3];
  Float_t z = 0.;
  Float_t Rp[2][3]={{0.,0.,-22.25},{0.,0.,-22.25}};
  Double_t weight;
  T->Branch("nmb",&nmb,"nmb/I");
  T->Branch("Q2",&Q2,"Q2/F");
  T->Branch("Xb",&Xb,"Xb/F");
  T->Branch("nu",&nu,"nu/F");
  T->Branch("q",q,"q[3]/F");
  T->Branch("Pe",Pe,"Pe[3]/F");
  T->Branch("Pe_size",&Pe_size,"Pe_size/F");
  T->Branch("theta_e",&theta_e,"theta_e/F");
  T->Branch("phi_e",&phi_e,"phi_e/F");
  T->Branch("Pp",Pp,"Pp[nmb][3]/F");
  T->Branch("Pp_size",Pp_size,"Pp_size[nmb]/F");
  T->Branch("pq_angle",pq_angle,"pq_angle[nmb]/F");
  T->Branch("Ep",Ep,"Ep[nmb]/F");
  T->Branch("theta_p",theta_p,"theta_p[nmb]/F");
  T->Branch("phi_p",phi_p,"phi_p[nmb]/F");
  T->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[nmb]/F");
  T->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[nmb]/F");
  T->Branch("Pmiss",Pmiss,"Pmiss[nmb][3]/F");
  T->Branch("Rp",Rp,"Rp[nmb][3]/F");
  T->Branch("weight",&weight,"weight/D");

  // Loop over all events
  const int nEvents = inTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %10000==0)
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);

      // Require a leading proton
      if (lead_type != pCode)
	continue;

      // Create vectors for the particles
      TVector3 ve(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vq = TVector3(0.,0.,eg2beam) - ve;
      TVector3 vlead(gen_pLead[0],gen_pLead[1],gen_pLead[2]);
      TVector3 vrec(gen_pRec[0],gen_pRec[1],gen_pRec[2]);
      TVector3 vmiss=vlead-vq;
      TVector3 vcm=vmiss+vrec;
      TVector3 vrel=0.5*(vmiss-vrec);

      double gen_pMiss_Mag = vmiss.Mag();
      double gen_pe_Mag = ve.Mag();
      double gen_QSq = 2. * eg2beam * gen_pe_Mag * (1. - ve.CosTheta());
      double gen_nu = eg2beam - ve.Mag();
      double gen_xB = gen_QSq/(2.*mN*gen_nu);
      double gen_q_Mag = vq.Mag();
      double gen_pLead_Mag = vlead.Mag();
      double gen_pRec_Mag = vrec.Mag();

      // Apply weight for detecting e, p      
      weight = gen_weight * eMap.accept(ve) * pMap.accept(vlead) * 1.E33; // put it in nb to make it macroscopic

      if (weight <= 0.)
	continue;

      // Do leading proton cuts
      if (vmiss.Mag() <0.3)
	continue;
      if (vmiss.Mag() >1.0)
	continue;
      if (gen_xB < 1.2)
	continue;
      if (vlead.Angle(vq)  > 25.*M_PI/180.)
	continue;
      if (gen_pLead_Mag/gen_q_Mag < 0.62)
	continue;
      if (gen_pLead_Mag/gen_q_Mag > 0.96)
	continue;
      if (sqrt(-gen_QSq + 4.*mN*gen_nu - 2.*sqrt(mN*mN + gen_pLead_Mag*gen_pLead_Mag)*(gen_nu + 2.*mN) + 5.*mN*mN + 2.*vq.Dot(vlead)) > 1.1)
	continue;
      if (!accept_electron(ve)) // Fiducial cut on electron
	continue;

      // Fill the acceptance histograms
      const double recoil_accept = pMap.accept(vrec);
      if (rec_type == pCode)
	{
	  h_rec_p_all->Fill(gen_pMiss_Mag,weight);
	  
	  // Test if the recoil was in the fiducial region and above threshold
	  if (accept_proton(vrec) && (vrec.Mag() > 0.35))
	    h_rec_p_acc->Fill(gen_pMiss_Mag,weight*recoil_accept);
	}

      // Decide which tree to write
      if (nmb == 1)
	{
	  // If the recoil is a proton, we better not have detected it.
	  if (rec_type == pCode)
	    weight *= (1. - recoil_accept);
	}
      else // nmb==2;
	{
	  if (rec_type == pCode)
	    weight *= recoil_accept;
	  else
	    weight = 0.; // We can't use recoil neutrons.
	}      

      // Load up tree
      Q2 = gen_QSq;
      Xb = gen_xB;
      nu = gen_nu;
      Pe_size = gen_pe_Mag;
      theta_e = ve.Theta() * 180./M_PI;
      phi_e = ve.Phi()*180./M_PI;
      if (phi_e < -30.) phi_e += 360;
      Pp_size[0] = gen_pLead_Mag;
      Pp_size[1] = gen_pRec_Mag;
      pq_angle[0] = vq.Angle(vlead)*180./M_PI;
      pq_angle[1] = vq.Angle(vrec)*180./M_PI;
      Ep[0] = sqrt(gen_pLead_Mag*gen_pLead_Mag + mN*mN);
      Ep[1] = sqrt(gen_pRec_Mag*gen_pRec_Mag + mN*mN);
      theta_p[0] = vlead.Theta()*180./M_PI;
      theta_p[0] = vrec.Theta()*180./M_PI;
      phi_p[0] = vlead.Phi()*180./M_PI;
      if (phi_p[0] < -30.) phi_p[0] += 360;
      phi_p[1] = vrec.Phi()*180./M_PI;
      if (phi_p[1] < -30.) phi_p[1] += 360;
      for (int i=0 ; i<3 ; i++)
	{
	  q[i] = vq[i];
	  Pe[i] = gen_pe[i];
	  Pp[0][i] = vlead[i];
	  Pp[1][i] = vrec[i];
	  Pmiss[0][i] = vlead[i] - vq[i];
	  Pmiss[1][i] = vrec[i] - vq[i];
	}

      Pmiss_q_angle[0] = (vlead - vq).Angle(vq) * 180./M_PI;
      Pmiss_q_angle[1] = (vrec - vq).Angle(vq) * 180./M_PI;
      
      Pmiss_size[0] = (vlead-vq).Mag();
      Pmiss_size[1] = (vrec - vq).Mag();

      T->Fill();
    }
  infile->Close();
  
  // Take the acceptance ratio
  //rec_p_rat->BayesDivide(h_rec_p_acc,h_rec_p_all);
  rec_p_rat->BayesDivide(h_rec_p_acc,h_rec_p_all,"cl=0.683 b(1,1) mode");

  outfile->cd();
  T->Write();
  rec_p_rat->Write();
  h_rec_p_all->Write();
  h_rec_p_acc->Write();

  outfile->Close();

  return 0;
}

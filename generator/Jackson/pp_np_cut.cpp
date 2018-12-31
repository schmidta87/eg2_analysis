#include <cstdlib>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include "AccMap.h"
#include "fiducials.h"
#include "Nuclear_Info.h"
#include "Cross_Sections.h"

using namespace std;

double sq(double x)
{
  return x*x;
}

int main(int argc, char ** argv)
{
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tsimulator /path/to/gen/file /path/to/file [# protons]\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  TFile * outfile = new TFile(argv[2],"RECREATE");
  double num_p = atof(argv[3]);
  
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
  TTree * T = new TTree("T","Cut Tree");
  Float_t Q2, Xb, Pe[3], Pe_size, theta_e, phi_e, Pp[2][3], Pp_size[2], pq_angle[2], Ep[2], theta_p[2], phi_p[2], nu, q[3];
  Float_t Pmiss_q_angle[2], Pmiss_size[2], Pmiss[2][3];
  Float_t z = 0.;
  Float_t Rp[2][3]={{0.,0.,-22.25},{0.,0.,-22.25}};
  Double_t weight;
  Double_t pRec_Mag;
  Double_t pLead_Mag;
  Int_t nmb = 2;
  
  T->Branch("Q2",&Q2,"Q2/F");
  T->Branch("Xb",&Xb,"Xb/F");
  T->Branch("nu",&nu,"nu/F");
  T->Branch("q",q,"q[3]/F");
  T->Branch("Pe",Pe,"Pe[3]/F");
  T->Branch("Pe_size",&Pe_size,"Pe_size/F");
  T->Branch("theta_e",&theta_e,"theta_e/F");
  T->Branch("phi_e",&phi_e,"phi_e/F");
  T->Branch("Pp",Pp,"Pp[2][3]/F");
  T->Branch("Pp_size",Pp_size,"Pp_size[2]/F");
  T->Branch("pq_angle",pq_angle,"pq_angle[2]/F");
  T->Branch("Ep",Ep,"Ep[2]/F");
  T->Branch("theta_p",theta_p,"theta_p[2]/F");
  T->Branch("phi_p",phi_p,"phi_p[2]/F");
  T->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[2]/F");
  T->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[2]/F");
  T->Branch("Pmiss",Pmiss,"Pmiss[2][3]/F");
  T->Branch("Rp",Rp,"Rp[2][3]/F");
  T->Branch("weight",&weight,"weight/D");
  T->Branch("pRec",&pRec_Mag,"pRec/D");
  T->Branch("pLead",&pLead_Mag,"pLead/D");

  TRandom3 myRand(0);
  double a = 0.03;
  double b = 0.06;
  Cross_Sections myCS; 
  
  // Loop over all events
  const int nEvents = inTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %100000==0)
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);
      
      // Require a recoil proton
      if (rec_type != pCode)
	continue;

      if ((lead_type == pCode) == (num_p == 1))
	continue;
      
      // Create vectors for the particles
      TVector3 vrec(gen_pRec[0],gen_pRec[1],gen_pRec[2]);
      TVector3 ve(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vq = TVector3(0.,0.,eg2beam) - ve;
      TVector3 vlead(gen_pLead[0],gen_pLead[1],gen_pLead[2]);
      
      // Fiducial cuts
      if (lead_type == nCode)
	vlead += vlead*myRand.Gaus(0,a+b*vlead.Mag());
      
      if (!accept_electron(ve))
	continue;
      if (!accept_proton(vlead))
	continue;
      if (!accept_proton(vrec))
	continue;
      if (!accept_neutron(vlead))
	continue;

      if (lead_type == pCode)
	vlead += vlead*myRand.Gaus(0,a+b*vlead.Mag());

      if (gen_weight <= 0.)
	continue;
      
      TVector3 vmiss=vlead-vq;
      TVector3 vcm=vmiss+vrec;
      TVector3 vrel=0.5*(vmiss-vrec);
      
      double gen_pMiss_Mag = vmiss.Mag();
      double gen_pe_Mag = ve.Mag();
      double gen_nu = eg2beam - ve.Mag();
      double gen_QSq = 2. * eg2beam * gen_pe_Mag * (1. - ve.CosTheta());
      double gen_xB = gen_QSq/(2.*mN*gen_nu);
      double gen_q_Mag = vq.Mag();
      double gen_pLead_Mag = vlead.Mag();
      double gen_pRec_Mag = vrec.Mag();
      double gen_ELead = sqrt(gen_pLead_Mag*gen_pLead_Mag+mN*mN);
      double m_miss = sqrt((gen_nu+2*mN-gen_ELead)*(gen_nu+2*mN-gen_ELead)-gen_pMiss_Mag*gen_pMiss_Mag);

      
      // Meytal's cuts
      if (gen_xB < 1.1)
	continue;
      if (0.62 > gen_pLead_Mag/gen_q_Mag)
	continue;
      if (gen_pLead_Mag/gen_q_Mag > 1.1)
	continue;
      if (vlead.Angle(vq)  > 25.*M_PI/180.)
	continue;
      if (m_miss > 1.175)
	continue;
      if (0.4 > gen_pMiss_Mag)
	continue;
      if (gen_pMiss_Mag > 1)
	continue;
      if (gen_pRec_Mag < 0.35)
	continue;

      // Apply weight for cross sections
      double gen_theta = ve.Theta();
      double tau = gen_QSq/(4*mN*mN);
      double epsilon = 1/(1.0+2.0*(1.+tau)*sq(tan(gen_theta/2)));
      
      if (lead_type == pCode)
	weight = gen_weight/(2*(epsilon/tau*sq(myCS.GEp(gen_QSq))+sq(myCS.GMp(gen_QSq))));
      else
	weight = gen_weight/(epsilon/tau*sq(myCS.GEn(gen_QSq))+sq(myCS.GMn(gen_QSq)));

      // Load up tree
      pRec_Mag = gen_pRec_Mag;
      pLead_Mag = gen_pLead_Mag;

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
  
  outfile->cd();
  T->Write();
  outfile->Close();

  return 0;
}


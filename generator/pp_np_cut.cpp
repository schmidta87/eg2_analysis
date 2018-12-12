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
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tsimulator /path/to/gen/file /path/to/pp/file /path/tp/np/file\n";
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
  Double_t weight;
  Double_t pRec_Mag;
  T->Branch("weight",&weight,"weight/D");
  T->Branch("pRec",&pRec_Mag,"pRec/D");
  
  // Loop over all events
  const int nEvents = inTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %10000==0)
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);

      // Require a recoil proton
      if (rec_type != pCode)
	continue;

      // Apply weight for detecting e, p      
      // weight = gen_weight * eMap.accept(ve) * pMap.accept(vlead) * 1.E33; // put it in nb to make it macroscopic

      weight = gen_weight;
      if (weight <= 0.)
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

      // Fiducial cuts
      if (!accept_electron(ve))
	continue;
      if (!accept_proton(vlead))
	continue;
      if (!accept_proton(vrec))
	continue;
      if (!accept_neutron(vlead))
	continue;

      // Load up tree
      pRec_Mag = gen_pRec_Mag;
      
      if ((lead_type == pCode) == (num_p == 1))
	continue;
      
      T->Fill();
    }
  
  infile->Close();
  
  outfile->cd();
  T->Write();
  outfile->Close();

  return 0;
}


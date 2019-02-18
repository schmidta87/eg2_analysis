#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include "AccMap.h"
#include "fiducials.h"
#include "Nuclear_Info.h"

using namespace std;

double eSmearing=0.003;
double pSmearing=0.01;
double Tp = 0.53;
double sig_Tp = 0.05;
double Tpp = 0.44;
double sig_Tpp = 0.04;

int main(int argc, char ** argv)
{
  if (argc < 6)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tsimulator /path/to/gen/file /path/to/farm/file /path/to/map/file /path/to/1p/file /path/to/2p/file\n\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  TFile * farmfile = new TFile(argv[2]);
  AccMap pMap(argv[3]);
  AccMap eMap(argv[3],"e");

  bool verbose = false;
  bool rand_flag = false;
  bool doSmearing=true;

  int c;
  while ((c=getopt (argc-4, &argv[4], "vre:p:O")) != -1)
    switch(c)
      {
      case 'v':
	verbose = true;
	break;
      case 'r':
	rand_flag = true;
	break;
      case 'e':
	eSmearing = atof(optarg);
	break;
      case 'p':
	pSmearing = atof(optarg);
	break;
      case 'O':
	doSmearing = false;
	break;
      case '?':
	return -1;
      default:
	abort();
      }

  // Input Tree
  TTree * inTree = (TTree*)infile->Get("T");
  Double_t gen_pe[3], gen_pLead[3], gen_weight;
  inTree->SetBranchAddress("weight",&gen_weight);
  inTree->SetBranchAddress("pe",gen_pe);
  inTree->SetBranchAddress("pLead",gen_pLead);

  TTree * farmTree = (TTree*)farmfile->Get("data");
  Int_t gPart, particles[12];
  farmTree->SetBranchAddress("gPart",&gPart);
  farmTree->SetBranchAddress("particle",particles);

  // Other set up
  TRandom3 myRand(0);

  if (rand_flag)
    {
      eSmearing = 0.0025 + myRand.Uniform()*0.001;
      pSmearing = 0.008 + myRand.Uniform()*0.004;
      Tp += myRand.Gaus(0,sig_Tp);
      Tpp += myRand.Gaus(0,sig_Tpp);
      if (verbose)
	{
	  cout << "Electron resolution selected as " << eSmearing*100 << "%.\n";
	  cout << "Proton resolution selected as " << pSmearing*100 << "%.\n";
	  cout << "Lead proton transparency selected as " << Tp << ".\n";
	  cout << "Two-proton transparency selected as " << Tpp << ".\n";
	}
    }

  TFile * outfile1p = new TFile(argv[4],"RECREATE");

  // Output Tree (1p)
  TTree * T1p = new TTree("T","Simulated Data Tree");
  Float_t Q2, Xb, Pe[3], Pe_size, theta_e, phi_e, Pp[2][3], Pp_size[2], pq_angle[2], Ep[2], theta_p[2], phi_p[2], nu, q[3];
  Float_t Pmiss_q_angle[2], Pmiss_size[2], Pmiss[2][3];
  Float_t z = 0.;
  Float_t Rp[2][3]={{0.,0.,-22.25},{0.,0.,-22.25}};
  Double_t weight;

  T1p->Branch("Q2",&Q2,"Q2/F");
  T1p->Branch("Xb",&Xb,"Xb/F");
  T1p->Branch("nu",&nu,"nu/F");
  T1p->Branch("q",q,"q[3]/F");
  T1p->Branch("Pe",Pe,"Pe[3]/F");
  T1p->Branch("Pe_size",&Pe_size,"Pe_size/F");
  T1p->Branch("theta_e",&theta_e,"theta_e/F");
  T1p->Branch("phi_e",&phi_e,"phi_e/F");
  T1p->Branch("Pp",Pp,"Pp[1][3]/F");
  T1p->Branch("Pp_size",Pp_size,"Pp_size[1]/F");
  T1p->Branch("pq_angle",pq_angle,"pq_angle[1]/F");
  T1p->Branch("Ep",Ep,"Ep[1]/F");
  T1p->Branch("theta_p",theta_p,"theta_p[1]/F");
  T1p->Branch("phi_p",phi_p,"phi_p[1]/F");
  T1p->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[1]/F");
  T1p->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[1]/F");
  T1p->Branch("Pmiss",Pmiss,"Pmiss[1][3]/F");
  T1p->Branch("Rp",Rp,"Rp[1][3]/F");
  T1p->Branch("weight",&weight,"weight/D");

  // Loop over all events (1p)
  if (verbose)
    cerr << "Looping over 1p events...\n";
  const int nEvents = farmTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %1000000==0 and verbose) 
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);
      farmTree->GetEvent(event);

      // Create vectors for the particles
      TVector3 ve(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vlead(gen_pLead[0],gen_pLead[1],gen_pLead[2]);
      
      // Smearing
      if (doSmearing)
	{
	  ve *= (1. + eSmearing * myRand.Gaus() * ve.Mag());
	  gen_pe[0] = ve.X();
	  gen_pe[1] = ve.Y();
	  gen_pe[2] = ve.Z();

	  vlead *= (1. + pSmearing * myRand.Gaus() * vlead.Mag());
	  gen_pLead[0] = vlead.X();
	  gen_pLead[1] = vlead.Y();
	  gen_pLead[2] = vlead.Z();
	}

      // Derived vectors
      TVector3 vq = TVector3(0.,0.,eg2beam) - ve;
      TVector3 vmiss=vlead-vq;

      double gen_pMiss_Mag = vmiss.Mag();
      double gen_pe_Mag = ve.Mag();
      double gen_QSq = 2. * eg2beam * gen_pe_Mag * (1. - ve.CosTheta());
      double gen_nu = eg2beam - ve.Mag();
      double gen_xB = gen_QSq/(2.*mN*gen_nu);
      double gen_q_Mag = vq.Mag();
      double gen_pLead_Mag = vlead.Mag();

      // Apply weight for detecting e, p      
      //weight = gen_weight * eMap.accept(ve) * pMap.accept(vlead) * 1.E33; // put it in nb to make it macroscopic
      weight = gen_weight * pMap.accept(vlead) * Tp; // put it in nb to make it macroscopic

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

      // Load up tree
      Q2 = gen_QSq;
      Xb = gen_xB;
      nu = gen_nu;
      Pe_size = gen_pe_Mag;
      theta_e = ve.Theta() * 180./M_PI;
      phi_e = ve.Phi()*180./M_PI;
      if (phi_e < -30.) phi_e += 360;
      Pp_size[0] = gen_pLead_Mag;
      pq_angle[0] = vq.Angle(vlead)*180./M_PI;
      Ep[0] = sqrt(gen_pLead_Mag*gen_pLead_Mag + mN*mN);
      theta_p[0] = vlead.Theta()*180./M_PI;
      phi_p[0] = vlead.Phi()*180./M_PI;
      if (phi_p[0] < -30.) phi_p[0] += 360;
      for (int i=0 ; i<3 ; i++)
	{
	  q[i] = vq[i];
	  Pe[i] = gen_pe[i];
	  Pp[0][i] = vlead[i];
	  Pmiss[0][i] = vlead[i] - vq[i];
	}

      Pmiss_q_angle[0] = (vlead - vq).Angle(vq) * 180./M_PI;
      
      Pmiss_size[0] = (vlead-vq).Mag();

      T1p->Fill();
    }

  outfile1p->cd();
  T1p->Write();
  outfile1p->Close();

  TFile * outfile2p = new TFile(argv[5],"RECREATE");

  // Output Tree (2p)
  TTree * T2p = new TTree("T","Simulated Data Tree");
  T2p->Branch("Q2",&Q2,"Q2/F");
  T2p->Branch("Xb",&Xb,"Xb/F");
  T2p->Branch("nu",&nu,"nu/F");
  T2p->Branch("q",q,"q[3]/F"); 
  T2p->Branch("Pe",Pe,"Pe[3]/F");
  T2p->Branch("Pe_size",&Pe_size,"Pe_size/F");
  T2p->Branch("theta_e",&theta_e,"theta_e/F");
  T2p->Branch("phi_e",&phi_e,"phi_e/F");
  T2p->Branch("Pp",Pp,"Pp[2][3]/F");
  T2p->Branch("Pp_size",Pp_size,"Pp_size[2]/F");
  T2p->Branch("pq_angle",pq_angle,"pq_angle[2]/F");
  T2p->Branch("Ep",Ep,"Ep[2]/F");
  T2p->Branch("theta_p",theta_p,"theta_p[2]/F");
  T2p->Branch("phi_p",phi_p,"phi_p[2]/F");
  T2p->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[2]/F");
  T2p->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[2]/F");
  T2p->Branch("Pmiss",Pmiss,"Pmiss[2][3]/F");
  T2p->Branch("Rp",Rp,"Rp[2][3]/F");
  T2p->Branch("weight",&weight,"weight/D");

  outfile2p->cd();
  T2p->Write();
  outfile2p->Close();

  infile->Close();
  farmfile->Close();

  return 0;
}
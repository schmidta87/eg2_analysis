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

AccMap * pMap;

// Hardcoded stuff
double eSmearing=0.003;
double pSmearing=0.01;
double Tp = 0.53;
double sig_Tp = 0.05;
double Tpp = 0.44;
double sig_Tpp = 0.04;

bool passes_lead_cuts(TVector3 ve, TVector3 vp);

int main(int argc, char ** argv)
{
  if (argc < 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tmistaken_recoils /path/to/gen/file /path/to/map/file /path/to/out/file\n\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  pMap = new AccMap(argv[2]);

  bool quiet = false;
  bool rand_flag = false;
  bool doSmearing = true;
  bool doTrans = true;
  bool doMaps = true;
  bool doFCuts = true;
  bool doSCuts = true;
  int c;
  while ((c=getopt (argc-3, &argv[3], "qre:p:OoMCS")) != -1)
    switch(c)
      {
      case 'q':
	quiet = true;
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
      case 'o':
	doTrans = false;
	break;
      case 'M':
	doMaps = false;
	break;
      case 'C':
	doFCuts = false;
	break;
      case 'S':
	doSCuts = false;
	break;
      case '?':
	return -1;
      default:
	abort();
      }

  // Histograms
  TH1D * h_pmiss_correct = new TH1D("pmiss_cor","Correct Leads;pmiss [GeV];Counts",28,0.3,1.0);
  h_pmiss_correct->Sumw2();
  TH1D * h_pmiss_mistake = new TH1D("pmiss_mis","Correct Leads;pmiss [GeV];Counts",28,0.3,1.0);
  h_pmiss_mistake->Sumw2();
  TH1D * h_pmiss_ratio;

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
  
  // Other set up
  TRandom3 myRand(0);

  if (rand_flag)
    {
      eSmearing = 0.0025 + myRand.Uniform()*0.001;
      pSmearing = 0.008 + myRand.Uniform()*0.004;
      Tp += myRand.Gaus(0,sig_Tp);
      Tpp += myRand.Gaus(0,sig_Tpp);
      if (!quiet)
	{
	  cout << "Electron resolution selected as " << eSmearing*100 << "%.\n";
	  cout << "Proton resolution selected as " << pSmearing*100 << "%.\n";
	  cout << "Lead proton transparency selected as " << Tp << ".\n";
	  cout << "Two-proton transparency selected as " << Tpp << ".\n";
	}
    }
  if (!doTrans){
    Tp=1;
    Tpp=1;
    if (!quiet){
      cout << "Transparency has been turned off.\n";
    }
  }

  TFile * outfile = new TFile(argv[3],"RECREATE");

  // Loop over all events (1p)
  const int nEvents = inTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if ((event %1000000==0) && (!quiet))
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);

      // Create vectors for the particles
      TVector3 ve(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vlead(gen_pLead[0],gen_pLead[1],gen_pLead[2]);
      TVector3 vrec(gen_pRec[0],gen_pRec[1],gen_pRec[2]);
      TVector3 vq = TVector3(0.,0.,eg2beam) - ve;

      // Do smearing
      if (doSmearing)
	{
	  ve *= (1. + eSmearing * myRand.Gaus() * ve.Mag());
	  vlead *= (1. + pSmearing * myRand.Gaus() * vlead.Mag());
	  vrec *= (1. + pSmearing * myRand.Gaus() * vrec.Mag());
	}
      
      // Test if each ep combo passes lead cuts
      if ((lead_type == pCode) && passes_lead_cuts(ve,vlead))
	h_pmiss_correct->Fill((vlead - vq).Mag(),1.E33 * gen_weight * pMap->accept(vlead));
      
      if ((rec_type == pCode) && passes_lead_cuts(ve,vrec))
	h_pmiss_mistake->Fill((vrec - vq).Mag(),1.E33 * gen_weight * pMap->accept(vrec));
    }

  infile->Close();
  outfile->cd();

  // Form the ratio histogram
  h_pmiss_ratio = (TH1D*) h_pmiss_mistake->Clone("pmiss_rat");
  h_pmiss_ratio->Divide(h_pmiss_correct);
  h_pmiss_ratio->SetTitle("Fraction of mistaken leads");

  // Write out
  h_pmiss_correct->Write();
  h_pmiss_mistake->Write();
  h_pmiss_ratio->Write();
  outfile->Close();

  return 0;
}


bool passes_lead_cuts(TVector3 ve, TVector3 vp)
{
  // Define q and q-related cuts  
  TVector3 q = TVector3(0.,0.,eg2beam) - ve;

  if (q.Angle(vp) > 25. *M_PI/180.)
    return false;

  double p_over_q = vp.Mag() / q.Mag();

  if (p_over_q < 0.62)
    return false;

  if (p_over_q > 0.96)
    return false;

  // Define QSq, xB
  double omega = eg2beam - ve.Mag();
  double QSq = q.Mag2() - omega*omega;
  double xB = QSq / (2.*mN*omega);
  
  if (xB < 1.2)
    return false;

  // Define pmiss
  TVector3 pmiss = vp - q;

  if (pmiss.Mag() < 0.3) 
    return false;

  if (pmiss.Mag() > 1.0)
    return false;

  // Missing mass cut
  if (sqrt(-QSq + 4.*mN*omega - 2.*sqrt(mN*mN + vp.Mag2())*(omega + 2.*mN) + 5.*mN*mN + 2.*q.Dot(vp)) > 1.1)
    return false;

  // Other cuts in make_hists
  if (vp.Mag() > 2.4)
    return false;

  if (!accept_electron(ve))
    return false;

  if (!accept_proton(vp))
    return false;

  // Apply an additional fiducial cut that the map acc must be > 80%                                                                   
  if ( pMap->accept(vp) < 0.8)
    return false;

  // This proton passes leading cuts!
  return true;
}

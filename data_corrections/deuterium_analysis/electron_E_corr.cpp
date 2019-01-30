#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"

using namespace std;

double expected_E(TVector3 ve, TVector3 vp);
double sq(double x){return x*x;};

int main(int argc, char ** argv)
{
  // Set up the input file
  TFile * inFile = new TFile("e1p_corr.root");
  TTree * inTree = (TTree*)inFile->Get("T");
  const int nEvents = inTree->GetEntries();
  Float_t electron_z, Pe[3], Pp[1][3], Rp[1][3], phi_e;
  inTree->SetBranchAddress("electron_z",&electron_z);
  inTree->SetBranchAddress("Pe",Pe);
  inTree->SetBranchAddress("phi_e",&phi_e);
  inTree->SetBranchAddress("Pp",Pp);
  inTree->SetBranchAddress("Rp",Rp);

  // Create histograms
  TFile * outFile = new TFile("output.root","RECREATE");
  TH2D *corr_by_theta[6], *corr_by_mom[6];
  for (int i=0 ; i<6 ; i++)
    {
      char name[100];
      char title[100];

      sprintf(name,"corr_by_theta_%d",i);
      sprintf(title,"Correction for sector %d;Theta [deg];E_meas - E_exp [GeV];Counts",i);
      corr_by_theta[i] = new TH2D(name,title,50,0.,50.,50,-0.05,0.05);

      sprintf(name,"corr_by_mom_%d",i);
      sprintf(title,"Correction for sector %d;Mom [GeV];E_meas - E_exp [GeV];Counts",i);
      corr_by_mom[i] = new TH2D(name,title,50,0.,5.,50,-0.05,0.05);
    }

  // Loop over events
  for (int event=0 ; event <nEvents ; event++)
    {
      if ((event%100000)==0)
        cerr << "Working on event " << event << " out of " << nEvents << "\n";
      inTree->GetEvent(event);

      if (fabs(Rp[0][2] - electron_z) > 1.5)
        continue;

      TVector3 ve(Pe[0],Pe[1],Pe[2]);
      TVector3 vp(Pp[0][0], Pp[0][1], Pp[0][2]);

      int sec = phi_e/60.;
      double delta_Ee = ve.Mag() - expected_E(ve,vp);

      corr_by_theta[sec]->Fill(ve.Theta()*180./M_PI,delta_Ee);
      corr_by_mom[sec]->Fill(ve.Mag(),delta_Ee);
    }

  inFile->Close();
  outFile->cd();
  for (int i=0 ; i<6 ; i++)
    {
      corr_by_mom[i]->Write();
      corr_by_theta[i]->Write();
    }
  outFile->Close();
  return 0;
}

const double mP = 0.9383;
const double mN = 0.9396;
const double mD = 1.8756;
const double eg2beam=5.0094;

double expected_E(TVector3 ve, TVector3 vp)
{
  double pp = vp.Mag();
  double Ep = sqrt(sq(pp) + sq(mP));
  double A = sq(mD) + sq(mP) - 2.*eg2beam*mD - 2.*eg2beam*pp + 2.*mD*pp + 2.*eg2beam*vp.Z();
  double B = -2.*eg2beam + 2.*mD + 2.*Ep + 2.*vp.Dot(ve.Unit()) + 2.*eg2beam*ve.CosTheta();

  return (sq(mN) - A)/B;
}

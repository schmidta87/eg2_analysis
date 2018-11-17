#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

#include "constants.h"
#include "fiducials.h"

using namespace std;

double constrainedGauss(double *x, double *p);
double specialLine(double *x, double *p);

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "/path/to/2p/file /path/to/out/file\n\n";
      return -1;
    }

  TFile * outfile = new TFile(argv[2],"RECREATE");

  // Initialize some histograms
  TH2D * hist_epp_cm_lon = new TH2D("cm_lon","epp events;p_miss [GeV];p_cm_lon [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_inp = new TH2D("cm_inp","epp events;p_miss [GeV];p_cm_inp [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_oop = new TH2D("cm_oop","epp events;p_miss [GeV];p_cm_oop [GeV];Counts",7,0.3,1.,24,-1.2,1.2);

  // Read in epp file
  cerr << "Loading 2p data file " << argv[2] << " ...\n";
  Float_t pPVec[2][3]; // Proton momenta
  Float_t qVec[3]; // q
  Float_t pMissVec[2][3];
  Float_t vertices[2][3]; // position 3 vectors of proton vertices
  TFile *fepp = new TFile(argv[1]);
  TTree *tepp = (TTree*)fepp->Get("T");
  tepp->SetBranchAddress("q",qVec);
  tepp->SetBranchAddress("Pmiss",pMissVec);
  tepp->SetBranchAddress("Rp",vertices);
  tepp->SetBranchAddress("Pp",pPVec);

  // See if there is a weight branch
  Double_t weight = 1.;
  TBranch * weight_branch = t1p->GetBranch("weight");
  if (weight_branch)
    {
      tepp->SetBranchAddress("weight",&weight);
    }

  for (int i=0 ; i<tepp->GetEntries() ; i++)
    {
      tepp->GetEvent(i);  

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 prec(pPVec[1][0],pPVec[1][1],pPVec[1][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);

      // Make Or's cut on leading proton
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
	continue;

      // Make Or's cut on recoil proton
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > min_prec)))
	continue;

      // Apply fiducial cuts, which Or has not yet made
      if (!accept_proton(prec))
	continue;

      double pmiss_mag = pmiss.Mag();

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = oop.Cross(lon).Unit();

      TVector3 pcm = prec + pmiss;
      double pcm_lon = pcm.Dot(lon);
      double pcm_inp = pcm.Dot(inp);
      double pcm_oop = pcm.Dot(oop);

      hist_epp_cm_lon->Fill(pmiss_mag,pcm_lon,weight);
      hist_epp_cm_inp->Fill(pmiss_mag,pcm_inp,weight);
      hist_epp_cm_oop->Fill(pmiss_mag,pcm_oop,weight);

    }
  fepp->Close();

  cerr << "Done reading input data.\n";

  // Fit the inp and oop histograms
  TH1D * histSumInpOop =  hist_epp_cm_inp->ProjectionY("histSumInpOop");
  TH1D * histTemp = hist_epp_cm_oop->ProjectionY("temp");
  TF1 * fInpOop = new TF1("fInpOop",constrainedGauss,-1.2,1.2,2);
  fInpOop->SetParameter(0,histSumInpOop->GetMaximum());
  fInpOop->SetParameter(1,0.1);
  histSumInpOop->Add(histTemp);
  histSumInpOop->Fit("fInpOop","Q");

  // Loop over slices
  double bins[5];
  double binErrs[5];
  double means_inp[5];
  double meanErrs_inp[5];
  double sigs_inp[5];
  double sigErrs_inp[5];
  double means_oop[5];
  double meanErrs_oop[5];
  double sigs_oop[5];
  double sigErrs_oop[5];
  double means_lon[5];
  double meanErrs_lon[5];
  double sigs_lon[5];
  double sigErrs_lon[5];


  TF1 * fGauss = new TF1("fGauss","gaus(0)",-1.2,1.2);
  int startBins[5]={1,2,3,4,5};
  int stopBins[5]={1,2,3,4,7};
  for (int i=0 ; i<5 ; i++)
    {
      bins[i] = 0.5*(hist_epp_cm_lon->GetXaxis()->GetBinCenter(startBins[i]) + 
			    hist_epp_cm_lon->GetXaxis()->GetBinCenter(stopBins[i]) );
      binErrs[i] =0.;
      char tempString[100];

      // First, the important longitudinal direction
      // Create the projection
      sprintf(tempString,"temp_%d",i);
      TH1D * histLon =  hist_epp_cm_lon->ProjectionY(tempString,startBins[i],stopBins[i]);

      // Set the initial guesses
      fGauss->SetParameter(0,histLon->GetMaximum());
      fGauss->SetParameter(1,histLon->GetBinCenter(histLon->GetMaximumBin()));
      fGauss->SetParameter(2,0.2);

      // Fit and extract results
      histLon->Fit("fGauss","Q");
      means_lon[i] = fGauss->GetParameter(1);
      meanErrs_lon[i] = fGauss->GetParError(1);
      sigs_lon[i] = fabs(fGauss->GetParameter(2));
      sigErrs_lon[i] = fGauss->GetParError(2);

      // Let's artificially blow up the errors if there aren't that many counts
      if (histLon->Integral()<5)
	{
	  meanErrs_lon[i]=10.*fabs(means_lon[i]);
	  sigErrs_lon[i]=10.*fabs(sigs_lon[i]);
	}

      // Next, the inp direction
      sprintf(tempString,"temp_inp_%d",i);
      TH1D * histInp =  hist_epp_cm_inp->ProjectionY(tempString,startBins[i],stopBins[i]);

      // Set the initial guesses
      fGauss->SetParameter(0,histInp->GetMaximum());
      fGauss->SetParameter(1,0.01);
      fGauss->SetParameter(2,0.2);

      // Fit and extract results
      histInp->Fit("fGauss","Q");
      means_inp[i] = fGauss->GetParameter(1);
      meanErrs_inp[i] = fGauss->GetParError(1);
      sigs_inp[i] = fabs(fGauss->GetParameter(2));
      sigErrs_inp[i] = fGauss->GetParError(2);

      // Next, the oop direction
      sprintf(tempString,"temp_oop_%d",i);
      TH1D * histOop =  hist_epp_cm_oop->ProjectionY(tempString,startBins[i],stopBins[i]);

      // Set the initial guesses
      fGauss->SetParameter(0,histOop->GetMaximum());
      fGauss->SetParameter(1,0.01);
      fGauss->SetParameter(2,0.2);

      // Fit and extract results
      histOop->Fit("fGauss","Q");
      means_oop[i] = fGauss->GetParameter(1);
      meanErrs_oop[i] = fGauss->GetParError(1);
      sigs_oop[i] = fabs(fGauss->GetParameter(2));
      sigErrs_oop[i] = fGauss->GetParError(2);
    }

  // Now do a fit to the longitudinal means and sigmas!
  TF1 * fLine = new TF1("fLine",specialLine,-1.2,1.2,2);
  TGraphErrors * meanGraph_lon = new TGraphErrors(5,bins,means_lon,binErrs,meanErrs_lon);
  TGraphErrors * sigGraph_lon = new TGraphErrors(5,bins,sigs_lon,binErrs,sigErrs_lon);
  TGraphErrors * meanGraph_inp = new TGraphErrors(5,bins,means_inp,binErrs,meanErrs_inp);
  TGraphErrors * sigGraph_inp = new TGraphErrors(5,bins,sigs_inp,binErrs,sigErrs_inp);
  TGraphErrors * meanGraph_oop = new TGraphErrors(5,bins,means_oop,binErrs,meanErrs_oop);
  TGraphErrors * sigGraph_oop = new TGraphErrors(5,bins,sigs_oop,binErrs,sigErrs_oop);
  fLine->SetParameter(0,0.1);
  fLine->SetParameter(1,0.1);
  meanGraph_lon->Fit("fLine","Q");
  double b1 = fLine->GetParameter(0);
  double b1_err = fLine->GetParError(0);
  double b2 = fLine->GetParameter(1);
  double b2_err = fLine->GetParError(1);
  fLine->SetParameter(0,0.1);
  fLine->SetParameter(1,0.1);
  sigGraph_lon->Fit("fLine","Q");
  double a1 = fLine->GetParameter(0);
  double a1_err = fLine->GetParError(0);
  double a2 = fLine->GetParameter(1);
  double a2_err = fLine->GetParError(1);

  // Print out fit results
  cout << "\n\na1 was fit to be " << a1 << " +/- " << a1_err << "\n";
  cout << "a2 was fit to be " << a2 << " +/- " << a2_err << "\n";
  cout << "b1 was fit to be " << b1 << " +/- " << b1_err << "\n";
  cout << "b2 was fit to be " << b2 << " +/- " << b2_err << "\n";
  cout << "SigPerp was fit to be " << fInpOop->GetParameter(1) << " +/- " << fInpOop->GetParError(1) << "\n\n\n";

  // Write out histograms
  outfile->cd();
  hist_epp_cm_lon->Write();
  hist_epp_cm_inp->Write();
  hist_epp_cm_oop->Write();
  histSumInpOop->Write();
  meanGraph_lon->Write("means_lon");
  sigGraph_lon->Write("sigmas_lon");
  meanGraph_inp->Write("means_inp");
  sigGraph_inp->Write("sigmas_inp");
  meanGraph_oop->Write("means_oop");
  sigGraph_oop->Write("sigmas_oop");
  outfile->Close();

  // Clean up
  return 0;
}

double constrainedGauss(double *x, double *p)
{
  double amp = p[0];
  double sig = p[1];
  double xOverSig = *x/sig;
  return amp*exp(-0.5 * xOverSig*xOverSig);
}

double specialLine(double *x, double *p)
{
  double a1 = p[0];
  double a2 = p[1];
  return a1*(*x - 0.6) + a2;
}

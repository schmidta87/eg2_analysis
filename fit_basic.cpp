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
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > 0.35)))
	continue;

      double pmiss_mag = pmiss.Mag();

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = lon.Cross(TVector3(0.,0.,1.)).Unit();
      TVector3 inp = lon.Cross(oop).Unit();

      TVector3 pcm = prec + pmiss;
      double pcm_lon = pcm.Dot(lon);
      double pcm_inp = pcm.Dot(inp);
      double pcm_oop = pcm.Dot(oop);

      hist_epp_cm_lon->Fill(pmiss_mag,pcm_lon);
      hist_epp_cm_inp->Fill(pmiss_mag,pcm_inp);
      hist_epp_cm_oop->Fill(pmiss_mag,pcm_oop);

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
  double means[5];
  double meanErrs[5];
  double sigs[5];
  double sigErrs[5];
  TF1 * fGauss = new TF1("fGauss","gaus(0)",-1.2,1.2);
  int startBins[5]={1,2,3,4,5};
  int stopBins[5]={1,2,3,4,7};
  for (int i=0 ; i<5 ; i++)
    {
      bins[i] = 0.5*(hist_epp_cm_lon->GetXaxis()->GetBinCenter(startBins[i]) + 
			    hist_epp_cm_lon->GetXaxis()->GetBinCenter(stopBins[i]) );
      binErrs[i] =0.;

      // Create the projection
      char tempString[100];
      sprintf(tempString,"temp_%d",i);
      TH1D * histLon =  hist_epp_cm_lon->ProjectionY(tempString,startBins[i],stopBins[i]);

      // Set the initial guesses
      fGauss->SetParameter(0,histLon->GetMaximum());
      fGauss->SetParameter(1,histLon->GetBinCenter(histLon->GetMaximumBin()));
      fGauss->SetParameter(2,0.2);

      // Fit and extract results
      histLon->Fit("fGauss","Q");
      means[i] = fGauss->GetParameter(1);
      meanErrs[i] = fGauss->GetParError(1);
      sigs[i] = fabs(fGauss->GetParameter(2));
      sigErrs[i] = fGauss->GetParError(2);

      // Let's artificially blow up the errors if there aren't that many counts
      if (histLon->Integral()<5)
	{
	  meanErrs[i]=10.*fabs(means[i]);
	  sigErrs[i]=10.*fabs(sigs[i]);
	}
    }

  // Now do a fit to the longitudinal means and sigmas!
  TF1 * fLine = new TF1("fLine",specialLine,-1.2,1.2,2);
  TGraphErrors * meanGraph = new TGraphErrors(5,bins,means,binErrs,meanErrs);
  TGraphErrors * sigGraph = new TGraphErrors(5,bins,sigs,binErrs,sigErrs);
  fLine->SetParameter(0,0.1);
  fLine->SetParameter(1,0.1);
  meanGraph->Fit("fLine","Q");
  double b1 = fLine->GetParameter(0);
  double b1_err = fLine->GetParError(0);
  double b2 = fLine->GetParameter(1);
  double b2_err = fLine->GetParError(1);
  fLine->SetParameter(0,0.1);
  fLine->SetParameter(1,0.1);
  sigGraph->Fit("fLine","Q");
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
  meanGraph->Write("means");
  sigGraph->Write("sigmas");
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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TVector3.h"

#include "constants.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use: "
	   << "\tlc_a2_test /path/to/generator/file /path/to/output/file [Q2 cut (GeV^2)]\n\n";
      return -1;
    }

  // Helper constants
  const double Ebeam=eg2beam;
  //const double Ebeam=10.;

  // Set up the input branches
  TFile * infile = new TFile(argv[1]);
  TTree * inT = (TTree*)infile->Get("genT");
  double pe[3], weight_nr, weight_lc;
  inT->SetBranchAddress("pe",pe);
  inT->SetBranchAddress("weight",&weight_nr);
  inT->SetBranchAddress("lcweight",&weight_lc);

  // Set up the outfile
  ofstream outfile(argv[2]);

  // xB histograms
  const double QSqcut=atof(argv[3]);
  TH1D * h_nr = new TH1D("nr","Non-relativistic;xB;Counts",40,1.,2.);
  h_nr->Sumw2();
  TH1D * h_lc = new TH1D("lc","Light-cone;xB;Counts",40,1.,2.);
  h_lc->Sumw2();

  // Event loop
  const int nEvents = inT->GetEntries();
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event % 10000==0)
	cerr << "Working on event " << event << " out of " << nEvents << "...\n";
      
      inT->GetEvent(event);

      // Assemble variables.
      TVector3 ve(pe[0],pe[1],pe[2]);
      TVector3 q = TVector3(0.,0.,Ebeam) - ve;
      double omega = Ebeam - ve.Mag();
      double QSq = q.Mag2() - omega*omega;
      double xB = QSq/(2.*mN*omega);

      //cerr << "xB: " << xB << "   QSq = " << QSq << "\n";

      // Cut on Q2
      if (QSq < QSqcut)
	continue;

      h_nr->Fill(xB,weight_nr*1.E33);
      h_lc->Fill(xB,weight_lc*1.E33);
    }
  
  for (int bin=1 ; bin<=h_nr->GetXaxis()->GetNbins() ; bin++)
    {
      outfile << h_nr->GetBinCenter(bin) << " "
      	      << h_nr->GetBinContent(bin) << " " << h_nr->GetBinError(bin) << " "
	      << h_lc->GetBinContent(bin) << " " << h_lc->GetBinError(bin) << "\n";
    }

  infile->Close();
  outfile.close();
  return 0;
}

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc != 5)
    {
      cerr << "Wrong number of arguments! Instead use\n"
	   << "\t/path/to/tree/file /path/to/output/text/file [# burn-in] [# to keep one of] \n\n";
      return -1;
    }

  // Get the parameters
  int burn_in = atoi(argv[3]);
  int keep_one = atoi(argv[4]);

  // Read in the file
  TFile * infile = new TFile(argv[1]);
  TTree * intree = (TTree*) infile->Get("mcmc");
  double a1, a2, b1, b2, sigPerp, logP;
  intree->SetBranchAddress("a1",&a1);
  intree->SetBranchAddress("a2",&a2);
  intree->SetBranchAddress("b1",&b1);
  intree->SetBranchAddress("b2",&b2);
  intree->SetBranchAddress("sigPerp",&sigPerp);
  intree->SetBranchAddress("logposterior",&logP);

  // Create the output file
  TFile * outfile = new TFile(argv[2],"RECREATE");
  TTree * outtree = new TTree("mcmc","Processed mcmc results");
  outtree->Branch("a1",&a1,"a1/D");
  outtree->Branch("a2",&a2,"a2/D");
  outtree->Branch("b1",&b1,"b1/D");
  outtree->Branch("b2",&b2,"b2/D");
  outtree->Branch("sigPerp",&sigPerp,"sigPerp/D");
  outtree->Branch("logposterior",&logP,"logP/D");

  // Loop over the tree
  for (int i=0 ; i<intree->GetEntries() ; i++)
    {
      intree->GetEvent(i);

      if (i<burn_in)
	continue;

      if (i%keep_one ==0)
	outtree->Fill();
    }

  outfile->cd();
  outtree->Write();
  outfile->Close();
  infile->Close();

  return 0;
}

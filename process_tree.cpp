#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 5)
    {
      cerr << "Wrong number of arguments! Instead use\n"
	   << "\t[# burn-in] [# to keep one of] /path/to/output/file List: /path/to/first/input /path/to/second/input ...\n\n";
      return -1;
    }

  // Get the parameters
  int burn_in = atoi(argv[1]);
  int keep_one = atoi(argv[2]);

  // Create the output file
  double a1, a2, b1, b2, sigPerp, logP;
  TFile * outfile = new TFile(argv[3],"RECREATE");
  TTree * outtree = new TTree("mcmc","Processed mcmc results");
  outtree->Branch("a1",&a1,"a1/D");
  outtree->Branch("a2",&a2,"a2/D");
  outtree->Branch("b1",&b1,"b1/D");
  outtree->Branch("b2",&b2,"b2/D");
  outtree->Branch("sigPerp",&sigPerp,"sigPerp/D");
  outtree->Branch("logposterior",&logP,"logP/D");

  // Loop over the input files
  for (int j=4 ; j<argc ; j++)
    {
      // Read in the file
      TFile * infile = new TFile(argv[j]);
      TTree * intree = (TTree*) infile->Get("mcmc");
      intree->SetBranchAddress("a1",&a1);
      intree->SetBranchAddress("a2",&a2);
      intree->SetBranchAddress("b1",&b1);
      intree->SetBranchAddress("b2",&b2);
      intree->SetBranchAddress("sigPerp",&sigPerp);
      intree->SetBranchAddress("logposterior",&logP);
      intree->SetBranchAddress("recoil_acc",&recoil_acc);

      // Loop over the tree
      for (int i=0 ; i<intree->GetEntries() ; i++)
	{
	  intree->GetEvent(i);
	  
	  if (i<burn_in)
	    continue;
	  
	  if (i%keep_one ==0)
	    outtree->Fill();
	}

      infile->Close();
    }
      
  outfile->cd();
  outtree->Write();
  outfile->Close();

  return 0;
}

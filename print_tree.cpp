#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments! Instead use\n"
	   << "\t/path/to/tree/file /path/to/output/text/file \n\n";
      return -1;
    }

  // Read in the file
  TFile * f = new TFile(argv[1]);
  TTree * t = (TTree*) f->Get("mcmc");
  double a1, a2, b1, b2, sigPerp, logP;
  t->SetBranchAddress("a1",&a1);
  t->SetBranchAddress("a2",&a2);
  t->SetBranchAddress("b1",&b1);
  t->SetBranchAddress("b2",&b2);
  t->SetBranchAddress("sigPerp",&sigPerp);
  t->SetBranchAddress("logposterior",&logP);

  // Create the output file
  ofstream outfile(argv[2]);

  // Loop over the tree
  for (int i=0 ; i<t->GetEntries() ; i++)
    {
      t->GetEvent(i);

      // Print out the entries
      outfile << a1 << " " << a2 << " " << b1 << " " << b2 << " " << sigPerp << " " << logP << "\n";
    }
  outfile.close();
  f->Close();

  return 0;
}

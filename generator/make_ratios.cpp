#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tmake_ratios /path/to/hist_file /path/to/outfile\n\n";
      return -1;
    }

  TFile * inf = new TFile(argv[1]);
  TH2D * h_ep = (TH2D*)inf->Get("ep_pmiss_E1");
  TH2D * h_epp = (TH2D*)inf->Get("epp_pmiss_E1");

  ofstream outf(argv[2]);
  outf << "# [pmiss (GeV/c)] [E1 (GeV)] [epp events] [ep events]\n";
  for (int binx=1 ; binx <= h_ep->GetXaxis()->GetNbins() ; binx++)
    {
      double pmiss = h_ep->GetXaxis()->GetBinCenter(binx);

      for (int biny=1 ; biny <= h_ep->GetYaxis()->GetNbins() ; biny++)
	{
	  double E1 = h_ep->GetYaxis()->GetBinCenter(biny);
	  outf << pmiss << " " << E1 << " " << h_epp->GetBinContent(binx,biny) << " " << h_ep->GetBinContent(binx,biny) << "\n";
	}
      outf << "\n";
    }
  outf.close();
  
  return 0;
}

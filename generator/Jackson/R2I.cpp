#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <ctype.h>
#include <fstream>
#include <string>

#include <sys/stat.h>
#include <sys/types.h>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1D.h"

using namespace std;

double M = 0.938;
int pCode = 2212;
int nCode = 2112;
int events_per_file = 2000;

int main(int argc, char ** argv)
{
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Instead try:\n\n ./R2I [path/to/input.root] [path/to/output/directory]\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  FILE *fp;
  std::string dir = argv[2];
     // Input Tree
  TTree * inTree = (TTree*)infile->Get("genT");
  Double_t pe[3], pLead[3], pRec[3], weight;
  Int_t lead_type, rec_type;
  inTree->SetBranchAddress("lead_type",&lead_type);
  inTree->SetBranchAddress("rec_type",&rec_type);
  inTree->SetBranchAddress("weight",&weight);
  inTree->SetBranchAddress("pe",pe);
  inTree->SetBranchAddress("pLead",pLead);
  inTree->SetBranchAddress("pRec",pRec);

  int nEvents = inTree->GetEntries();

  int nFiles = nEvents/events_per_file;
  if (nFiles > 2000)
    nFiles = 2000;

  for (int file=0 ; file < nFiles ; file ++)
    {
      char filestr[5];
      snprintf(filestr,5,"%d",file);
      std::string s = dir + "/" + filestr;
      const char *dirname = s.c_str();
      mkdir(dirname, 0755);
      
      std::string s2 = s + "/GenSource.inp";
      const char *filename = s2.c_str();
      cout << filename << "\n";
      fp = fopen(filename,"w");
      for (int event = 0 ; event < events_per_file ; event++)
	{
	  inTree->GetEvent(event + file * events_per_file);
	  int Zlead = 0;
	  int Zrec = 0;
	  if (lead_type == pCode)
	    Zlead = 1;
	  if (rec_type == pCode)
	    Zrec = 1;

	  fprintf(fp,"%4d %4d %8f %20f %20f %20f %20f %20f %20f %10d",
		 1,Zlead,M,0.0,0.0,0.0,pLead[0],pLead[1],pLead[2],event+1);
	  fprintf(fp,"\n");
	  fprintf(fp,"%4d %4d %8f %20f %20f %20f %20f %20f %20f %10d",
		 1,Zrec,M,0.0,0.0,0.0,pRec[0],pRec[1],pRec[2],event+1);
	  fprintf(fp,"\n");

	}
      fclose(fp);
    }

  return 0;
}

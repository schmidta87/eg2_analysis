#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TFile.h"
#include "TChain.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

using namespace std;

int main(int argc, char ** argv){

  if (argc != 2)
    {
      cerr << "Wrong number of arguments. Instead try\n\t"
	   <<"Emiss path/to/output/text/file";
      return -1;
    }
  //create output text file
  ofstream textfile;
  textfile.open(argv[1]);

  //work with root files in specific directory
  const char *dirname="/home/awild/nobackup/Data/eg2Project/poolEmiss";
  const char *ext=".root";
  //  const char *type="_11_"; fname.Contains(type,TString::kExact)
  TSystemDirectory dir(dirname, dirname); 
  TList *files = dir.GetListOfFiles(); 

  //make histograms for means and stds
  int nbins = 100;
  TH2D * pRec_eMiss_means = new TH2D("pRec_eMiss_means","pRec_eMiss;pRec;eMiss_means;Counts",nbins,0.35,0.9,nbins,0,0.5);
  TH2D * pRec_eMiss_stds = new TH2D("pRec_eMiss_stds","pRec_eMiss;pRec;eMiss_stds;Counts",nbins,0.35,0.9,nbins,0,0.125);
  //  pRec_eMiss_means->Sumw2();
  //  pRec_eMiss_stds->Sumw2();

  //This is the loop over files
  if (files) { 
    TSystemFile *file; 
    TString fname;
    TString fullname;
    TIter next(files); 
    while ((file=(TSystemFile*)next())) { 
      fname = file->GetName(); 
      fullname = "/home/awild/nobackup/Data/eg2Project/poolEmiss/" + fname;
      if (!file->IsDirectory() && fname.EndsWith(ext)) {

	
	//This is the loop over pRec Projections
	TFile * histFile = new TFile(fullname); 
	TH2D * inputHist = (TH2D*)histFile->Get("pRec_eMiss");
	for(int i = 0; i < inputHist->GetXaxis()->GetNbins(); i++){
	  TH1D * HistProj = inputHist->ProjectionY("pRec_Proj",i,(i+1));
	  //cout<<inputHist->GetXaxis()->GetBinCenter(i+1)<<" "<<HistProj->GetMean()<<" "<<HistProj->GetStdDev()<<endl;
	  pRec_eMiss_means->Fill(inputHist->GetXaxis()->GetBinCenter(i+1),HistProj->GetMean());
	  pRec_eMiss_stds->Fill(inputHist->GetXaxis()->GetBinCenter(i+1),HistProj->GetStdDev());
	}
	
	cout << fname.Data() << endl; 
      }     
    }
  }

  textfile<<"[Column 1: pRec Value (Means)] [Column 2: Mean Emiss Value (Means)] [Column 3: Emiss Stat Error (Means)] [Column 4: pRec Value (STD)] [Column 5: Emiss STD Value (STD)] [Column 3: Emiss Stat Error (STD)]"<<endl;
  for(int j = 0; j < nbins; j++){
    TH1D * meanHistProj = pRec_eMiss_means->ProjectionY("pRec_Proj_means",j,(j+1));
    TH1D * stdHistProj = pRec_eMiss_stds->ProjectionY("pRec_Proj_stds",j,(j+1));
    textfile << pRec_eMiss_means->GetXaxis()->GetBinCenter(j+1) << " " << meanHistProj->GetMean() << " " << meanHistProj->GetStdDev() << " " << pRec_eMiss_stds->GetXaxis()->GetBinCenter(j+1) << " " << stdHistProj->GetMean() << " " << stdHistProj->GetStdDev() << endl;
    
}

  textfile.close();
  return 0;
}

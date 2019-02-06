#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
using namespace std;

int main(int argc, char ** argv){
  if( argc != 3){
    cerr<<"Wrong number of arguments. Instead try:\n\t"
        << "H2T /path/to/input/root/file /path/to/output/text/file\n";
    return -1;
  }
  //Get data files                                                                                                                                                                                                 
  TFile * Data = new TFile(argv[1]);
  ofstream file;
  file.open(argv[2]);
  TList* List = (TList*)Data->GetListOfKeys();
  string keys[List->GetEntries()];
  cout<<"Keys Retrieved"<<endl;
  const string endchar = "goatsMaToats";
  bool isTH1D = false;
  bool isTH2D = false;
  for(int j = -1; j < 60; j++){
    cout<<"Entry number:"<<j<<endl;
    int k = 0;
    for(int i = 0; i < (List->GetEntries()); i++){
      isTH1D = false;//Data->Get(List->At(i)->GetName())->InheritsFrom(TH1D::Class());
      isTH2D = Data->Get(List->At(i)->GetName())->InheritsFrom(TH2D::Class());

      if(j>-1){


	if(isTH1D){     	  
	  TH1D * WHist = (TH1D*)Data->Get(List->At(i)->GetName());	  
	  if(j < WHist->GetXaxis()->GetNbins()){
	    file << WHist->GetBinCenter(j) << " " << WHist->GetBinContent(j) << " " << WHist->GetBinError(j) << " ";
	  }
	  else file <<"N N N ";
	}


	else if(isTH2D){
	  TH2D * SHist = (TH2D*)Data->Get(List->At(i)->GetName());
	  if(j < SHist->GetXaxis()->GetNbins()){
	    TH1D * HistProj = SHist->ProjectionY("epsilonProj",j,(j+1));
	    file << SHist->GetXaxis()->GetBinCenter(j) << " " << HistProj->GetMean() << " " << HistProj->GetStdDev() << " ";
	  }
	  else file <<"N N N ";
	  

	}
      }
      else{
	if(isTH1D || isTH2D){     
	  file <<"[Column "<<3*k+1<<": "<<List->At(i)->GetName()<<" xAxis]"<<"[Column "<<3*k+2<<": "<<List->At(i)->GetName()<<" yAxis]"<<"[Column "<<3*k+3<<": "<<List->At(i)->GetName()<<" yError]";
	  k++;

	  if(isTH1D) file<<" TH1D"<<endl;
	  else if(isTH2D) file<<" TH2D"<<endl;
	}
      }
    }
    file << endl;
  }
  cout << "finished \n";
  return 0;
}

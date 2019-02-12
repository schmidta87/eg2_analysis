#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

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
  bool isTGAE = false;
  
  //Go over all possible x bins
  for(int j = -1; j < 300; j++){

    //For Titles
    int k = 1;
    
    //Go over each TObject in List
    for(int i = 0; i < (List->GetEntries()); i++){
      isTH1D = Data->Get(List->At(i)->GetName())->InheritsFrom(TH1D::Class());
      isTH2D = Data->Get(List->At(i)->GetName())->InheritsFrom(TH2D::Class());
      isTGAE = Data->Get(List->At(i)->GetName())->InheritsFrom(TGraphAsymmErrors::Class());
      
      //For TH1D
      if(isTH1D){
	TH1D * WHist = (TH1D*)Data->Get(List->At(i)->GetName());	  
	string a = List->At(i)->GetName();
	int N = WHist->GetXaxis()->GetNbins();
	if(j<0){
	  file <<"[Column "<<k<<": "<<a<<" xAxis] [Column "<<(k+1)<<": "<<a<<" yAxis] [Column "<<(k+2)<<": "<<a<<" ySTD] \t TH1D"<<endl;
	  k=k+3;
	}
	else if((j>=0)&&(j<N)){
	  file << WHist->GetBinCenter(j) << " " << WHist->GetBinContent(j) << " " << WHist->GetBinError(j) << " ";
	}
	else if(j>=N) file <<"N N N ";
      }
      
    
      //For TH2D
      else if(isTH2D){
	TH2D * SHist = (TH2D*)Data->Get(List->At(i)->GetName());
	string a = List->At(i)->GetName();
	int N = SHist->GetXaxis()->GetNbins();
	if(j<0){
	  file <<"[Column "<<k<<": "<<a<<" xAxis] [Column "<<(k+1)<<": "<<a<<" yAxis] [Column "<<(k+2)<<": "<<a<<" ySTD] [Column "<<(k+3)<<": "<<a<<" yError] \t TH2D"<<endl;
	  k=k+4;      
	}
	else if((j>=0)&&(j<N)){
	  TH1D * HistProj = SHist->ProjectionY("epsilonProj",j,(j+1));
	  double n = HistProj->GetEntries();
	  double std = HistProj->GetStdDev();
	  file << SHist->GetXaxis()->GetBinCenter(j) << " " << HistProj->GetMean() << " " <<std<<" "<<(std/sqrt(n)) << " ";
	}
	else if(j>=N) file <<"N N N N ";
      }
      
      //For TGAE
      else if(isTGAE){
	TGraphAsymmErrors * myTGAE = (TGraphAsymmErrors*)Data->Get(List->At(i)->GetName());
	string a = List->At(i)->GetName();
	int N =myTGAE->GetN();
	if(j<0){
	  file <<"[Column "<<k<<": "<<a<<" xAxis] [Column "<<(k+1)<<": "<<a<<" xErrorLow] [Column "<<(k+2)<<": "<<a<<" xErrorHigh]"
	       <<"[Column "<<(k+3)<<": "<<a<<" yAxis] [Column "<<(k+4)<<": "<<a<<" yErrorLow] [Column "<<(k+5)<<": "<<a<<" yErrorHigh] \t TGAE"<<endl;
	  k=k+6;      
	}
	if((j>=0)&&(j<N)){
	  Double_t x,y;
	  myTGAE->GetPoint(j,x,y);
	  file << x << " " << myTGAE->GetErrorXlow(j) << " " << myTGAE->GetErrorXhigh(j) << " "
	       << y << " " << myTGAE->GetErrorYlow(j) << " " << myTGAE->GetErrorYhigh(j) << " ";
	}
	else if(j>=N) file <<"N N N N N N ";
      }
      
    }
    file << endl;
  }
  cout << "finished \n";
  return 0;
}

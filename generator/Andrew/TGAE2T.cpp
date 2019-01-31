#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
using namespace std;

int main(int argc, char ** argv){
  if( argc != 3){
    cerr<<"Wrong number of arguments. Instead try:\n\t"
        << "TGAE2T /path/to/input/root/file /path/to/output/text/file\n";
    return -1;
  }
  //Get data files                    
  TFile * Data = new TFile(argv[1]);
  ofstream file;
  file.open(argv[2]);
  const double coarse_bin_edges[11]={0.3,0.375,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.875,1.0};
  TH1D * px = (TH1D*)Data->Get("ep_Pm");
  TH1D * ppx = (TH1D*)Data->Get("epp_Pm");
  TH1D * newpx = (TH1D*) px->Rebin(10,"newpx",coarse_bin_edges);
  TH1D * newppx = (TH1D*) ppx->Rebin(10,"newppx",coarse_bin_edges);
  TGraphAsymmErrors * WGraph = new TGraphAsymmErrors;
  WGraph->BayesDivide(newppx,newpx);
  for(int j = -1; j < 60; j++){
    if(j>-1){
        if(j < WGraph->GetN()){
          Double_t x,y;
          WGraph->GetPoint(j,x,y);
          file << x << " " << WGraph->GetErrorXlow(j) << " " << WGraph->GetErrorXhigh(j) << " "
               << y << " " << WGraph->GetErrorYlow(j) << " " << WGraph->GetErrorYhigh(j) << " \n";
        }
        else{
          file <<"N N N N N N \n";
        }
      }
      else{
        file <<"[Column "<<1<<": pp_to_p xAxis]"
             <<"[Column "<<2<<": pp_to_p xLowerError]"
             <<"[Column "<<3<<": pp_to_p xHigherError]"
             <<"[Column "<<4<<": pp_to_p yAxis]"
             <<"[Column "<<5<<": pp_to_p yLowerError]"
             <<"[Column "<<6<<": pp_to_p yHigherError]"<<endl;

    }
    file << endl;
  }
  cout << "finished \n";
  return 0;
}

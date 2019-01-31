#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
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
  for(int j = -1; j < 60; j++){
     cout<<"Entry number:"<<j<<endl;
    for(int i = 1; i < (List->GetEntries()); i++){
      if(j>-1){

        TH1D * WHist = (TH1D*)Data->Get(List->At(i)->GetName());
        if(j < WHist->GetXaxis()->GetNbins()){
          file << WHist->GetBinCenter(j) << " " << WHist->GetBinContent(j) << " " << WHist->GetBinError(j) << " ";
        }
        else{
          file <<"N N N ";
        }
      }
      else{
        file <<"[Column "<<3*(i-1)+1<<": "<<List->At(i)->GetName()<<" xAxis]"<<"[Column "<<3*(i-1)+2<<": "<<List->At(i)->GetName()<<" yAxis]"<<"[Column "<<3*(i-1)+3<<": "<<List->At(i)->GetName()<<" yError]"<<endl;
      }
    }
    file << endl;
  }
  cout << "finished \n";
  return 0;
}

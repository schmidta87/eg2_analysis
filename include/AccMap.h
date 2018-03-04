#ifndef __ACCMAP_H__
#define __ACCMAP_H__

#include "TVector3.h"

class TH3D;
class TFile;

const int nBinsMom=100;
const int nBinsCos=200;
const int nBinsPhi=360;

class AccMap
{
 public:
  AccMap(const char * filename);
  ~AccMap();
  double accept(TVector3 p);
  double recoil_accept(TVector3 p);

 private:
  TFile * mapfile;
  TH3D * genHist;
  TH3D * accHist;
  double accept_map(TVector3 p);
  double accept_fake(TVector3 p);
  
};

#endif

#include "AccMap.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TH3.h"

#include "fiducials.h"
#include "constants.h"

AccMap::AccMap(const char * filename)
{
  mapfile = new TFile(filename);
  genHist = (TH3D*) mapfile->Get("solid_p_gen");
  accHist = (TH3D*) mapfile->Get("solid_p_acc");
}

AccMap::AccMap(const char * filename, const char * particle)
{
  mapfile = new TFile(filename);
  char temp[100];
  sprintf(temp,"solid_%s_gen",particle);
  genHist = (TH3D*) mapfile->Get(temp);
  sprintf(temp,"solid_%s_acc",particle);
  accHist = (TH3D*) mapfile->Get(temp);
}

AccMap::~AccMap()
{
  mapfile->Close();
}

double AccMap::recoil_accept(TVector3 p)
{
  if (accept_proton(p))
    if (p.Mag() > min_prec)
      return accept(p);
  return 0.;
}

double AccMap::accept(TVector3 p)
{
  // return 1.;
  //return accept_fake(p);
   return accept_map(p);
};

double AccMap::accept_fake(TVector3 p)
{
  double phiDeg = p.Phi()*180./M_PI;
  if (phiDeg < -30.)
    phiDeg += 360.;
  
  int sector = int(phiDeg/60.);
  double san_phiDeg = phiDeg - 60.*sector;
  
  if (fabs(san_phiDeg)<15.)
    return 1.;
  else
    return 0.;
}

double AccMap::accept_map(TVector3 p)
{
  double mom = p.Mag();
  double cosTheta = p.CosTheta();
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.)
    phi += 360.;

  int momBin = genHist->GetXaxis()->FindBin(mom);
  int cosBin = genHist->GetYaxis()->FindBin(cosTheta);
  int phiBin = genHist->GetZaxis()->FindBin(phi);

  // WE HAVE TO SANITIZE MOMENTUM BECAUSE THE MAPS ARE SLIGHTLY INCOMPLETE
  // EVENTUALLY WE'LL WANT TO REMOVE THIS!
  if (momBin > 70) momBin=70;
  
  double gen=0.;
  double acc=0.;
  int exp_size=-1;

  while ((gen <= 0.)&&(exp_size<3))
    {
      exp_size++;
      gen = acc = 0.;
      for (int mB = momBin-exp_size; (mB<=momBin+exp_size)&&(mB<=nBinsMom) ; mB++)
	{
	  for (int cB = cosBin-exp_size; (cB<=cosBin+exp_size)&&(cB<=nBinsCos) ; cB++)
	    {
	      for (int pB = phiBin-exp_size; pB<=phiBin+exp_size ; pB++)
		{
		  int actualPB = pB;
		  if (actualPB < 1) actualPB+=nBinsPhi;
		  else if (actualPB > nBinsPhi) actualPB-=nBinsPhi;
		  
		  gen += genHist->GetBinContent(mB,cB,actualPB);
		  acc += accHist->GetBinContent(mB,cB,actualPB);
		}
	    }
	}
    }

  if (gen>0.)
    return acc/gen;
  else
    return 1.;
}

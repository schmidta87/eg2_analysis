#ifndef __CROSS_SECTIONS_H__
#define __CROSS_SECTIONS_H__

#include "TVector3.h"

class Cross_Sections
{
 public:
  Cross_Sections();
  double sigmaCCn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n);
  double sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigmaCC2(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  
 private:
  double Gdipole(double QSq);
  double sq(double x);
  double dot4(double x0, TVector3 x, double y0, TVector3 y);
  
};

#endif

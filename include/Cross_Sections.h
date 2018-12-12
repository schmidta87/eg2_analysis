#ifndef __CROSS_SECTIONS_H__
#define __CROSS_SECTIONS_H__

#include "TVector3.h"

const double mu_p=2.79;
const double mu_n=-1.91;

enum ffModel {dipole, kelly};

class Cross_Sections
{
 public:
  Cross_Sections();
  double sigmaCCn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n);
  double sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigmaCC2(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double GEp(double QSq);
  double GEn(double QSq);
  double GMp(double QSq);
  double GMn(double QSq);

 private:
  ffModel myModel;
  static double Gdipole(double QSq);
  static double Gkelly(double QSq,double a1, double b1, double b2, double b3);
  
};

#endif

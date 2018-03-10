#include <cmath>
#include "TVector3.h"

#include <iostream>

#include "helpers.h"

const Double_t kThetaPar0PiPlus[6] = {  7.00823   ,  5.5        ,  7.06596   ,  6.32763   ,  5.5       ,  5.5      };
const Double_t kThetaPar1PiPlus[6] = {  0.207249  ,  0.1        ,  0.127764  ,  0.1       ,  0.211012  ,  0.281549 };
const Double_t kThetaPar2PiPlus[6] = {  0.169287  ,  0.506354   , -0.0663754 ,  0.221727  ,  0.640963  ,  0.358452 };
const Double_t kThetaPar3PiPlus[6] = {  0.1       ,  0.1        ,  0.100003  ,  0.1       ,  0.1       ,  0.1      };
const Double_t kThetaPar4PiPlus[6] = {  0.1       ,  3.30779    ,  4.499     ,  5.30981   ,  3.20347   ,  0.776161 };
const Double_t kThetaPar5PiPlus[6] = { -0.1       , -0.651811   , -3.1793    , -3.3461    , -1.10808   , -0.462045 };

const Double_t kFidPar0Low0PiPlus[6] = {  25.      ,  25.      ,  25.     ,  25.       ,  25.       ,  25.       };
const Double_t kFidPar0Low1PiPlus[6] = { -12.      , -12.      , -12.     , -12        , -12        , -12        };
const Double_t kFidPar0Low2PiPlus[6] = {   1.64476 ,   1.51915 ,   1.1095 ,   0.977829 ,   0.955366 ,   0.969146 };
const Double_t kFidPar0Low3PiPlus[6] = {   4.4     ,   4.4     ,   4.4    ,   4.4      ,   4.4      ,   4.4      };

const Double_t kFidPar1Low0PiPlus[6] = {  4.        ,  4.   ,  2.78427 ,  3.58539 ,  3.32277   ,  4.      };
const Double_t kFidPar1Low1PiPlus[6] = {  2.        ,  2.   ,  2.      ,  1.38233 ,  0.0410601 ,  2.      };
const Double_t kFidPar1Low2PiPlus[6] = { -0.978469  , -2.   , -1.73543 , -2.      , -0.953828  , -2.      };
const Double_t kFidPar1Low3PiPlus[6] = {  0.5       ,  0.5  ,  0.5     ,  0.5     ,  0.5       ,  1.08576 };

const Double_t kFidPar0High0PiPlus[6] = {  25.      , 24.8096  , 24.8758  ,  25.       , 25.       , 25.      };
const Double_t kFidPar0High1PiPlus[6] = { -11.9735  , -8.      , -8.      , -12.       , -8.52574  , -8.      };
const Double_t kFidPar0High2PiPlus[6] = {  0.803484 ,  0.85143 ,  1.01249 ,   0.910994 ,  0.682825 ,  0.88846 };
const Double_t kFidPar0High3PiPlus[6] = {  4.40024  ,  4.8     ,  4.8     ,   4.4      ,  4.79866  ,  4.8     };

const Double_t kFidPar1High0PiPlus[6] = {  2.53606  ,  2.65468  ,  3.17084 ,  2.47156 ,  2.42349  ,  2.64394 };
const Double_t kFidPar1High1PiPlus[6] = {  0.442034 ,  0.201149 ,  1.27519 ,  1.76076 ,  1.25399  ,  0.15892 };
const Double_t kFidPar1High2PiPlus[6] = { -2.       , -0.179631 , -2.      , -1.89436 , -2.       , -2.      };
const Double_t kFidPar1High3PiPlus[6] = {  1.02806  ,  1.6      ,  0.5     ,  1.03961 ,  0.815707 ,  1.31013 };

double a(double mom, double p0, double p1, double p2, double p3)
{
  return p0 + p1*exp(p2*(mom-p3));
}

double b(double mom, double p0, double p1, double p2, double p3)
{
  return p0 + p1*mom*exp(p2*sq(mom-p3));
}

double theta_min(double mom, double p0, double p1, double p2, double p3, double p4, double p5)
{
  return p0 + p1/sq(mom) + p2*mom + p3/mom + p4*exp(mom*p5);
}

double deltaPhi(double theta, double a, double b, double thetaMin)
{
  return a*(1. - 1./((theta-thetaMin)/b + 1.));
}

bool accept_proton(TVector3 p)
{
  double mom = p.Mag();
  double theta = p.Theta() * 180./M_PI;
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) phi+= 360.;

  int sector = (phi+30.)/60.;

  double minTheta = theta_min(mom,kThetaPar0PiPlus[sector],
			      kThetaPar1PiPlus[sector],
			      kThetaPar2PiPlus[sector],
			      kThetaPar3PiPlus[sector],
			      kThetaPar4PiPlus[sector],
			      kThetaPar5PiPlus[sector]);

  // Double check theta
  if (theta < minTheta)
    return false;

  double phiCentral = 60.*sector;

  double aLow = a(mom,kFidPar0Low0PiPlus[sector],
		  kFidPar0Low1PiPlus[sector],
		  kFidPar0Low2PiPlus[sector],
		  kFidPar0Low3PiPlus[sector]);

  double aHigh = a(mom,kFidPar0High0PiPlus[sector],
		  kFidPar0High1PiPlus[sector],
		  kFidPar0High2PiPlus[sector],
		  kFidPar0High3PiPlus[sector]);

  double bLow = b(mom,kFidPar1Low0PiPlus[sector],
		  kFidPar1Low1PiPlus[sector],
		  kFidPar1Low2PiPlus[sector],
		  kFidPar1Low3PiPlus[sector]);

  double bHigh = b(mom,kFidPar1High0PiPlus[sector],
		  kFidPar1High1PiPlus[sector],
		  kFidPar1High2PiPlus[sector],
		  kFidPar1High3PiPlus[sector]);

  double deltaPhiLow = deltaPhi(theta, aLow, bLow, minTheta);
  double deltaPhiHigh = deltaPhi(theta, aHigh, bHigh, minTheta);

  //std::cout << mom << " " << theta << " " << phi << " " << sector << " " << aLow << " " << aHigh << "\n";

  if ((phi < phiCentral + deltaPhiHigh) && (phi > phiCentral - deltaPhiLow))
    return true;

  return false;
}

#include <cmath>
#include "TVector3.h"

#include <iostream>

#include "helpers.h"

// ############## PROTON FIDUCIAL CUT PARAMETERS ################
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


// ############## ELECTRON FIDUCIAL CUT PARAMETERS (starting with e) #################
// sector, side
const double eFidAp0[6][2]={{25.9392,25.142},{26.3123,26.546},{25.2953,26.4484},{28.,23.},{28.,23.2826},{28.,23.8565}};
const double eFidAp1[6][2]={{-12.,-12.},{-12.,-12.},{-12.,-12.},{-12.,-12.},{-10.802,-12.},{-11.9703,-12.}};
const double eFidAp2[6][2]={{0.453685,0.357039},{0.40864,0.425442},{0.486011,0.519846},{0.366217,0.562825},{0.203741,1.15862},{0.283665,1.4323}};
const double eFidAp3[6][2]={{4.4,4.4},{4.4,4.4},{4.4,4.4},{4.4,4.4},{4.40456,4.4},{4.4497,4.4}};

const double eFidBp0[6][2]={{2.53279,2.25185},{2.9569,3.27384},{1.98429,1.5},{2.29784,1.65084},{2.79198,1.92076},{1.84305,3.18777}};
const double eFidBp1[6][2]={{1.19464,0.527297},{2.,2.},{1.914,1.81489},{0.351108,1.91044},{0.534785,0.90012},{0.883945,0.539951}};
const double eFidBp2[6][2]={{-2.,-0.1},{-2.,-2.},{-0.5088,-0.273122},{-0.1,-0.574873},{-2.,-0.112026},{-0.1,-2.}};
const double eFidBp3[6][2]={{1.25872,1.55534},{0.74865,0.519786},{0.721026,0.5},{1.6,0.5},{0.720147,0.692653},{1.6,1.02089}};

const double eFidThetaMinp0[6]={15.,13.,13.,13.,14.6669,14.0658};
const double eFidThetaMinp1[6]={4.,1.99082,3.82167,1.50112,1.79097,2.20492};
const double eFidThetaMinp2[6]={-0.828672,-0.159112,-0.0647212,-1.3,-0.355746,-1.3};
const double eFidThetaMinp3[6]={2.51825,0.100001,0.100008,0.1,0.1,0.1};
const double eFidThetaMinp4[6]={2.15567,13.6076,9.26858,12.3802,14.1265,8.90341};
const double eFidThetaMinp5[6]={-0.1,-0.554978,-0.448555,-0.200086,-0.777645,-0.1};
//##############################################################3


bool accept_electron(TVector3 p)
{
  double mom = p.Mag();
  double theta = p.Theta() * 180./M_PI;
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) phi+= 360.;

  int sector = (phi+30.)/60.;

  double minTheta = theta_min(mom,eFidThetaMinp0[sector],
			      eFidThetaMinp1[sector],
			      eFidThetaMinp2[sector],
			      eFidThetaMinp3[sector],
			      eFidThetaMinp4[sector],
			      eFidThetaMinp5[sector]);

  // Double check theta
  if (theta < minTheta)
    return false;

  double phiCentral = 60.*sector;

  double aLow = a(mom,eFidAp0[sector][0],
		  eFidAp1[sector][0],
		  eFidAp2[sector][0],
		  eFidAp3[sector][0]);

  double aHigh = a(mom,eFidAp0[sector][1],
		  eFidAp1[sector][1],
		  eFidAp2[sector][1],
		  eFidAp3[sector][1]);

  double bLow = a(mom,eFidBp0[sector][0],
		  eFidBp1[sector][0],
		  eFidBp2[sector][0],
		  eFidBp3[sector][0]);

  double bHigh = a(mom,eFidBp0[sector][1],
		  eFidBp1[sector][1],
		  eFidBp2[sector][1],
		  eFidBp3[sector][1]);

  double deltaPhiLow = deltaPhi(theta, aLow, bLow, minTheta);
  double deltaPhiHigh = deltaPhi(theta, aHigh, bHigh, minTheta);

  //std::cout << mom << " " << theta << " " << phi << " " << sector << " " << aLow << " " << aHigh << "\n";

  if ((phi < phiCentral + deltaPhiHigh) && (phi > phiCentral - deltaPhiLow))
    return true;

  return false;
}

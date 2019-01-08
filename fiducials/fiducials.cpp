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

// For FidThetaMin calculation for electron
const Double_t kThetaPar0[6] = { 15        , 15        ,  15       , 15        ,  13       ,  13        };
const Double_t kThetaPar1[6] = { -0.425145 , -1.02217  , -0.7837   , -1.47798  ,   3.47361 ,   3.5714   };
const Double_t kThetaPar2[6] = { -0.666294 , -0.616567 , -0.673602 , -0.647113 ,  -0.34459 ,  -0.398458 };
const Double_t kThetaPar3[6] = {  5.73077  ,  5.51799  ,  8.05224  ,  7.74737  ,   8.45226 ,   9.54265  };
const Double_t kThetaPar4[6] = { 10.4976   , 14.0557   , 15.2178   , 16.7291   , -63.4556  , -22.649    };
const Double_t kThetaPar5[6] = { -1.13254  , -1.16189  , -2.08386  , -1.79939  ,  -3.3791  ,  -1.89746  };

// For parameter 0 of the FidPhiMin calculation for electron
const Double_t kFidPar0Low0[6] = {  25      ,  25        ,  25       ,  24.6345  ,  23.4731  ,  24.8599  };
const Double_t kFidPar0Low1[6] = { -12      , -12        , -12       , -12       , -12       , -12       };
const Double_t kFidPar0Low2[6] = {   0.5605 ,   0.714261 ,  0.616788 ,   0.62982 ,   1.84236 ,   1.00513 };
const Double_t kFidPar0Low3[6] = {  4.4     ,  4.4       ,  4.4      ,   4.4     ,   4.4     ,   4.4     };

// For parameter 1 of the FidPhiMin calculation for electron
const Double_t kFidPar1Low0[6] = {  2.1945   ,  4        ,  3.3352  ,  2.22769   ,  1.63143   ,  3.19807  };
const Double_t kFidPar1Low1[6] = {  1.51417  ,  1.56882  ,  2       ,  2         ,  1.90179   ,  0.173168 };
const Double_t kFidPar1Low2[6] = { -0.354081 , -2        , -2       , -0.760895  , -0.213751  , -0.1      };
const Double_t kFidPar1Low3[6] = {  0.5      ,  0.5      ,  1.01681 ,  1.31808   ,  0.786844  ,  1.6      };

// For parameter 0 of the FidPhiMax calculation for electron
const Double_t kFidPar0High0[6] = { 25       ,  25        ,  25        ,  25        ,  23.7067  ,  25       };
const Double_t kFidPar0High1[6] = { -8       , -10.3277   , -12        , -11.3361   , -12       , -11.4641  };
const Double_t kFidPar0High2[6] = {  0.479446 ,  0.380908 ,   0.675835 ,   0.636018 ,   2.92146 ,   0.55553 };
const Double_t kFidPar0High3[6] = {  4.8      ,  4.79964  ,   4.4      ,   4.4815   ,   4.4     ,   4.41327 };

// For parameter 1 of the FidPhiMax calculation for electron
const Double_t kFidPar1High0[6] = {  3.57349 ,  3.02279  ,  2.02102 ,  3.1948   ,  3.0934   ,  2.48828 };
const Double_t kFidPar1High1[6] = {  2       ,  0.966175 ,  2       ,  0.192701 ,  0.821726 ,  2       };
const Double_t kFidPar1High2[6] = { -2       , -2        , -1.70021 , -1.27578  , -0.233492 , -2       };
const Double_t kFidPar1High3[6] = {  0.5     ,  0.527823 ,  0.68655 ,  1.6      ,  1.6      ,  0.70261 };


bool accept_electron(TVector3 p)
{
  double mom = p.Mag();
  double theta = p.Theta() * 180./M_PI;
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) phi+= 360.;

  int sector = (phi+30.)/60.;

  double minTheta = theta_min(mom,kThetaPar0[sector],
			      kThetaPar1[sector],
			      kThetaPar2[sector],
			      kThetaPar3[sector],
			      kThetaPar4[sector],
			      kThetaPar5[sector]);

  // Double check theta
  if (theta < minTheta)
    return false;

  double phiCentral = 60.*sector;


  double aLow = a(mom,kFidPar0Low0[sector],
		  kFidPar0Low1[sector],
		  kFidPar0Low2[sector],
		  kFidPar0Low3[sector]);

  double aHigh = a(mom,kFidPar0High0[sector],
		  kFidPar0High1[sector],
		  kFidPar0High2[sector],
		  kFidPar0High3[sector]);

  double bLow = b(mom,kFidPar1Low0[sector],
		  kFidPar1Low1[sector],
		  kFidPar1Low2[sector],
		  kFidPar1Low3[sector]);

  double bHigh = b(mom,kFidPar1High0[sector],
		  kFidPar1High1[sector],
		  kFidPar1High2[sector],
		  kFidPar1High3[sector]);
  /*
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
  */

  double deltaPhiLow = deltaPhi(theta, aLow, bLow, minTheta);
  double deltaPhiHigh = deltaPhi(theta, aHigh, bHigh, minTheta);

  //std::cout << mom << " " << theta << " " << phi << " " << sector << " " << aLow << " " << aHigh << "\n";

  if ((phi < phiCentral + deltaPhiHigh) && (phi > phiCentral - deltaPhiLow))
    return true;

  return false;
}

// ############## NEUTRON FIDUCIAL CUT FOR ECAL #################

bool accept_neutron(TVector3 pm){
        double vtx_z = -22.25;
	double dist = 10.0;
        double r, pm_phi, x_rot, y_rot, angle, pm_theta, third_theta, d_intercept;

        const double ec_theta =  65.*M_PI/180.; // degrees
        const double z0       =563.1; // cm

        // Using law of sines to determine distance from vertex to EC
        pm_theta    = pm.Theta();
        third_theta = (M_PI - ec_theta - pm_theta);

        r = (z0-vtx_z)*sin(ec_theta)/sin(third_theta);

        double Xmiss = r*sin(pm.Theta())*cos(pm.Phi());
        double Ymiss = r*sin(pm.Theta())*sin(pm.Phi());

        pm_phi = 180./M_PI*pm.Phi();
        if (pm_phi < -30.) pm_phi += 360.;
        int sec = (int)(pm_phi+30)/60;
        if (sec>5) sec = 5;
        if (sec<0) sec = 0;

        angle = M_PI/180.*60.*(sec);

        x_rot = Xmiss*cos(angle) + Ymiss*sin(angle);
        y_rot = Xmiss*sin(angle) - Ymiss*cos(angle);

        d_intercept = dist/cos(atan(1.73));

        return((x_rot<390.-dist)&&(x_rot>1.73*fabs(y_rot)+55.+d_intercept));
}


// ############## NEUTRON FIDUCIAL CUT FOR TOF #################

bool accept_neutron_tof(TVector3 p)
{
  // This is our fiducial limit, 10 cm from the edge of each bar.
  const double distFromEdge=10.;

  // angle in degrees of the bar centers
  double thetaBars[57]={8.6,10.3,11.9,13.6,15.3,17,18.8,20.5,22.3,24,25.8,27.6,29.3,31.1,32.8,34.5,36.2,37.9,39.6,
			41.2,42.8,44.4,45.9,47.4,49.6,51.9,54.3,56.8,59.4,62,64.7,67.4,70.2,72.9,75.7,78.2,80.8,83.5,
			86.3,89.3,92.2,95.3,98.4,101.6,104.8,108.8,112,114.9,118.1,121.4,124.7,128,131.4,134.2,136.5,138.8,141}; 
  
  // width (length?) of the bars in cm
  double width[57]={32.3,48.1,64,79.8,95.7,106.6,122.4,138.3,154.1,170,185.8,201.7,217.6,233.4,249.3,265.1,281,296.8,312.7,
		    328.5,344.4,360.2,376.1,371.3,378.2,385,391.9,398.7,405.6,412.5,419.3,426.2,433,439.9,445.1,439.3,433.6,427.8,
		    422,416.3,410.5,404.8,399,393.3,387.5,380.1,363.5,347,330.4,313.9,297.3,280.8,264.2,246.8,235.4,224,212.7};

  // Distance to the target center in cm
  double distance[57]={513,509,503,500,498,496,495,494,493,493,493,494,495,496,498,501,504,507,511,
		       515,519,524,529,514,504,495,487,480,473,468,463,460,458,457,457,446,437,428,
		       421,414,409,405,402,400,399,402,395,389,384,381,378,377,378,379,380,383,385};
  
  
  double thetaDeg = p.Theta()*180./M_PI;
  if (thetaDeg<9 || thetaDeg>140)
    return false;
  
  // Sanitize phi to [-30,330]
  double phiDeg = p.Phi()*180./M_PI;
  if (phiDeg <-30.) phiDeg += 360.;

  // Get the value of phi within the relevant sector
  int sector = phiDeg/60;
  double PhiOneSector = phiDeg - 60.*sector;
  
  // Figure out the nearest bars
  int dsBar, usBar;
  for (dsBar=0, usBar=1; usBar < 57 ; dsBar++, usBar++)
    {
      if ((thetaBars[dsBar] < thetaDeg) && (thetaBars[usBar] > thetaDeg))
	break;
    }

  // Make a linear interpolation between the nearest bars
  double dsPhiMax = 180./M_PI*atan((width[dsBar] - distFromEdge)/(sin(thetaBars[dsBar]*M_PI/180.)*distance[dsBar]));
  double usPhiMax = 180./M_PI*atan((width[usBar] - distFromEdge)/(sin(thetaBars[usBar]*M_PI/180.)*distance[usBar]));
  double phiMax = dsPhiMax + (usPhiMax-dsPhiMax) * (thetaDeg - thetaBars[dsBar])/(thetaBars[usBar] - thetaBars[dsBar]);

  // Test if our neutron is within the phi bound
  if (fabs(PhiOneSector) < phiMax)
    return true;

  return false;
}

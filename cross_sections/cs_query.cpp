#include "Cross_Sections.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "TVector3.h"

using namespace std;

int main(int argc, char **argv)
{
  if (argc!=3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "  cs_query [csMethod: onshell, cc1, cc2] [ffModel: dipole, kelly]\n";
      return -1;
    }

  csMethod csMeth;
  ffModel ffMod;

  // Establish the Cross Section Method
  if (strcmp(argv[1],"onshell")==0)
    csMeth=onshell;
  else if (strcmp(argv[1],"cc1")==0)
    csMeth=cc1;
  else if (strcmp(argv[1],"cc2")==0)
    csMeth=cc2;
  else
    {
      cerr << "Invalid cs method!\n";
      return -2;
    }

  // Establish the Form Factor Model
  if (strcmp(argv[2],"dipole")==0)
    ffMod=dipole;
  else if (strcmp(argv[2],"kelly")==0)
    ffMod=kelly;
  else
    {
      cerr << "Invalid ff model!\n";
      return -3;
    }

  Cross_Sections myCS(csMeth,ffMod);  

  cerr << "Enter lines in the format: \n"
       << "[Ebeam GeV] [mom_e GeV] [theta_e deg] [phi_e deg] [mom_p GeV] [theta_p (deg)] [phi_p (deg)]\n";

  double Ebeam, mom_e, thetaDeg_e, phiDeg_e, mom_p, thetaDeg_p, phiDeg_p;
  while( cin >> Ebeam )
    {
      cin >> mom_e >> thetaDeg_e >> phiDeg_e >> mom_p >> thetaDeg_p >> phiDeg_p;

      TVector3 pe(mom_e * sin(thetaDeg_e * M_PI/180.) * cos(phiDeg_e*M_PI/180.),
		  mom_e * sin(thetaDeg_e * M_PI/180.) * sin(phiDeg_e*M_PI/180.),
		  mom_e * cos(thetaDeg_e * M_PI/180.));

      TVector3 pp(mom_p * sin(thetaDeg_p * M_PI/180.) * cos(phiDeg_p*M_PI/180.),
		  mom_p * sin(thetaDeg_p * M_PI/180.) * sin(phiDeg_p*M_PI/180.),
		  mom_p * cos(thetaDeg_p * M_PI/180.));

      // Calculate a few more things for Ronen
      TVector3 q = TVector3(0.,0.,Ebeam) - pe;
      double phiAngleDeg = q.Cross(pe).Angle(q.Cross(pp))*180./M_PI;
      double gammaDeg = q.Angle(pp) * 180./M_PI;

      cout << Ebeam << " " << mom_e << " " << thetaDeg_e << " " << phiDeg_e << " " << mom_p << " " <<thetaDeg_p << " " << phiDeg_p << " "
	   << " " << phiAngleDeg << " " << gammaDeg << " " << myCS.sigma_eN(Ebeam,pe,pp,true) << "\n";

    }
  return 0;
}

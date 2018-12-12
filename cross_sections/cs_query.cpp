#include "Cross_Sections.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "TVector3.h"

using namespace std;

int main(int argc, char **argv)
{
  Cross_Sections myCS;

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
		  

      cout << Ebeam << " " << mom_e << " " << thetaDeg_e << " " << phiDeg_e << " " << mom_p << " " <<thetaDeg_p << " " << phiDeg_p << " "
	   << myCS.sigmaCC1(Ebeam,pe,pp,true) << "\n";

    }
  return 0;
}

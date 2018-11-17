#include "fiducials.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "TVector3.h"
#include "TRandom3.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 2)
    {
      cerr << "Call fid_test using: \n"
	   << "\tfid_test /path/to/out/file\n\n";
      return -1;
    }

  ofstream outfile(argv[1]);
  TRandom3 myRand(0);

  for (int i=0 ; i < 10000 ; i++)
    {
      double theta = myRand.Rndm()*M_PI;
      double phi = myRand.Rndm()*2.*M_PI - M_PI/6.;
      double mom = 2.;

      TVector3 p(mom*sin(theta)*cos(phi), mom*sin(theta)*sin(phi), mom*cos(theta));
      int p_fid = accept_proton(p)? 1: 0;
      int e_fid = accept_electron(p)? 1: 0;

      outfile << theta*180./M_PI << " " << phi*180./M_PI << " " << e_fid << " " << p_fid << "\n";
    }

  outfile.close();

  return 0;
}

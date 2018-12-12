#include "AccMap.h"

#include <iostream>
#include <cmath>

#include "TVector3.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Call map_test using: \n"
	   << "\tmap_test /path/to/map/file [particle]\n\n";
      return -1;
    }

  AccMap myMap(argv[1],argv[2]);

  double mom, theta, phi;
  while (true)
    {
      cin >> mom >> theta >> phi;
      
      if (! cin.good())
        break;
      
      double x = mom*sin(theta*M_PI/180.)*cos(phi*M_PI/180.);
      double y = mom*sin(theta*M_PI/180.)*sin(phi*M_PI/180.);
      double z = mom*cos(theta*M_PI/180.);

      TVector3 p(x,y,z);
      cout << mom << " " << theta << " " << phi << " " << myMap.accept(p) << "\n";
    }

  return 0;
}

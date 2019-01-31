#ifndef __FIDUCIALS_H__
#define __FIDUCIALS_H__

#include "TVector3.h"

bool accept_proton(TVector3 p);
bool accept_electron(TVector3 p);
bool accept_neutron(TVector3 pm);
bool accept_neutron_tof(TVector3 p);

#endif

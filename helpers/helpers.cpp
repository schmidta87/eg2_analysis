#include "helpers.h"

double sq(double x)
{
  return x*x;
}

int clas_sector(double phi_deg)
{
  while (phi_deg < -30.)
    phi_deg+=360.;

  while (phi_deg > 330.)
    phi_deg-=360.;

  return (phi_deg+30)/6;
}

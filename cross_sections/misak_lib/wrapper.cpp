#include <cmath>
#include "misak.h"

extern"C"{
  void sigma_en_bound_(int * icon, int * nuc_id, float * ei, float *q2, float * p_ini, float * th_i, float * phi_i,
		       float * sigma_en_df1, float * sigma_en_st, float * sigma_en_free, float * sigma_en_lc,
		       float * eta, float * ifl);
}

// ************************************************************
// icon  = 0 - initialization, 1 - calculation
// nuc_id = 1 - proton -1 - neutron
// ei = initial electron energy GeV
// q2 = Q2
// p_ini = inital momenutm of the nucleon GeV/c
// th_i = angle of initial nucleon vs q - in radians
// phi_i = azimuthal angle of initial nucleon in the Ref frame with (qe) = (zx)
// sigma_en_df1  - de-forest CC1
// sigma_en_St   - for stational nucleon with Ei and theta_e
// sigma_en_free - for free nucleon with the same initial momentum
// sigma_en_lc   - light-cone approximation
// eta           -   virtuality parameter Eq.(30) Phys. Rev C08, 035202, 2018
//                   https://arxiv.org/abs/1805.04639
// ifl           - 0 - ok, 1 - kinematically forbidden
// ************************************************************

void misak_init()
{
  // Prepare arguments
  int icon = 0;
  int nuc_id = 1;
  float ei = 5.;
  float q2 = 2.;
  float p_ini = 0.4;
  float th_i = 0.2;
  float phi_i = 0.1;
  float sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc;
  float eta, ifl;

  sigma_en_bound_(&icon, &nuc_id, &ei, &q2, &p_ini, &th_i, &phi_i,
		  &sigma_en_df1, &sigma_en_st, &sigma_en_free, &sigma_en_lc,
		  &eta, &ifl);

}

void call_misak(double Ebeam, TVector3 k, TVector3 p, bool isProton,
		float &sigma_en_df1, float &sigma_en_st, float &sigma_en_free, float & sigma_en_lc)
{
  // Prep results
  sigma_en_df1=0.;
  sigma_en_st=0.;
  sigma_en_free=0.;
  sigma_en_lc=0.;

  // Helpers
  TVector3 vq = TVector3(0,0,Ebeam) - k;
  double nu = Ebeam - k.Mag();
  TVector3 vmiss = p - vq;

  // Misak's reference frame
  TVector3 z_hat = vq.Unit();
  TVector3 y_hat = vq.Cross(k).Unit();
  TVector3 x_hat = y_hat.Cross(z_hat);

  // Prepare arguments
  int icon = 1;
  int nuc_id = isProton ? 1 : -1;
  float ei = Ebeam;
  float q2 = vq.Mag2() - nu*nu;
  float p_ini = vmiss.Mag();
  float th_i = vmiss.Angle(vq);
  float phi_i = atan2(vmiss.Dot(y_hat),vmiss.Dot(x_hat));
  float eta, ifl;

  sigma_en_bound_(&icon, &nuc_id, &ei, &q2, &p_ini, &th_i, &phi_i,
		  &sigma_en_df1, &sigma_en_st, &sigma_en_free, &sigma_en_lc,
		  &eta, &ifl);
}

double sigma_misak_df1(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  float sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc;
  call_misak(Ebeam, k, p, isProton, sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc);
  return sigma_en_df1;
}

double sigma_misak_st(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  float sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc;
  call_misak(Ebeam, k, p, isProton, sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc);
  return sigma_en_st;
}

double sigma_misak_free(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  float sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc;
  call_misak(Ebeam, k, p, isProton, sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc);
  return sigma_en_free;
}

double sigma_misak_lc(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  float sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc;
  call_misak(Ebeam, k, p, isProton, sigma_en_df1, sigma_en_st, sigma_en_free, sigma_en_lc);
  return sigma_en_lc;
}

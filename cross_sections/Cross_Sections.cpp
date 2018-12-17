#include "Cross_Sections.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TVector3.h"

#include "constants.h"
#include "helpers.h"

Cross_Sections::Cross_Sections()
{
  myModel=kelly;
  myMethod=cc1;
}

Cross_Sections::Cross_Sections(csMethod thisMeth, ffModel thisMod)
{
  myModel=thisMod;
  myMethod=thisMeth;
}

Cross_Sections::~Cross_Sections(){}

double Cross_Sections::sigma_eN(double Ebeam,TVector3 k, TVector3 p, bool isProton)
{
  switch (myMethod)
    {
    case onshell:
      return sigma_onShell_by_Etheta(Ebeam,k,isProton);
    case cc1:
      return sigmaCC1(Ebeam,k,p,isProton);
    case cc2:
      return sigmaCC2(Ebeam,k,p,isProton);
    default:
      {
	std::cerr << "Invalid cross section method! Double check and fix!\n";
	exit(-1);
      }
    }
  return 0;
}

// This is the opposite from our usual (1,-1,-1,-1) convention
double dot4(double x0, TVector3 x, double y0, TVector3 y)
{
  return ((x*y)-(x0*y0));
}

double Cross_Sections::sigmaCCn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n)
{
  TVector3 q = TVector3(0.,0.,Ebeam) - k;
  TVector3 pM = p-q;
  
  double omega = Ebeam - k.Mag();
  double QSq = q.Mag2() - sq(omega);
  double E = sqrt(p.Mag2() + sq(mN));
  double Ebar = sqrt(pM.Mag2() + sq(mN));
  double omegabar = E-Ebar;
  double QSqbar = q.Mag2() - sq(omegabar);

  // Calculate form factors
  double GE = (isProton)? GEp(QSq) : GEn(QSq);
  double GM = (isProton)? GMp(QSq) : GMn(QSq);
  double F1 = 0.5 * (GE + QSq*GM/(4.*sq(mN)));
  double kF2 = (GM - GE)/(1. + QSq/(4.*sq(mN)));

  double wC;
  double wT;
  double wS;
  double wI;
  
  if (n==1)
    {
      wC = (sq(E+Ebar)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2)) - q.Mag2()*sq(F1 + kF2))/(4.*E*Ebar);
      wT = QSqbar*sq(F1 + kF2)/(2.*Ebar*E);
      wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
      wI = -p.Mag()*sin(p.Angle(q))*(Ebar + E)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
    }
  else if (n==2)
    {  
      double pbarp = dot4(Ebar,pM,E,p);
      double pbarq = dot4(Ebar,pM,omega,q);
      double pq = dot4(E,p,omega,q);
      double qbarq = dot4(omegabar,q,omega,q);
      double sumq = dot4((Ebar+E),(pM+p),omega,q);

      wC = ((E*Ebar + 0.5 * (pbarp + sq(mN))) * sq(F1)
	    - 0.5 * q.Mag2() * F1 * kF2
	    - ((pbarq*E + pq*Ebar)*omega
		      - Ebar * E * QSq
		      + pbarq * pq
		      - 0.5 * (pbarp - sq(mN)) * q.Mag2())
	    * sq(kF2)/(4*sq(mN)))
	/(E*Ebar);
      wT = (-(pbarp + sq(mN)) * sq(F1)
	    + qbarq * F1 * kF2
	    + (2*pbarq*pq - (pbarp - sq(mN))*QSq)
	    * sq(kF2)/(4*sq(mN)))
	/(Ebar*E);
      wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSq * sq(kF2) / (4*sq(mN)) ) / (E*Ebar);
      wI = p.Mag() * sin(p.Angle(q)) * (-(Ebar + E) * sq(F1)
					+ (sumq * omega - (Ebar + E) * QSq)
					* sq(kF2)/(4*sq(mN)))
	/(E*Ebar);
    }
  else
    {
      std::cerr << "Invalid cross section designation. Check and fix. Exiting\n\n\n";
    }
      
  double sigmaMott = cmSqGeVSq * 4. * sq(alpha) * k.Mag2() * sq(cos(k.Theta()/2.)) / sq(QSq);

  double cosPhi = cos(q.Cross(k).Angle( q.Cross(p) ));
  return sigmaMott * ( sq(QSq)/sq(q.Mag2()) * wC +
                       (QSq/(2.*q.Mag2()) + sq(tan(k.Theta()/2.))) * wT +
                       QSq/q.Mag2() * sqrt(QSq/q.Mag2() + sq(tan(k.Theta()/2.))) * wI * cosPhi +
                       (QSq/q.Mag2() * sq(cosPhi) + sq(tan(k.Theta()/2.))) * wS
                       );
}

double Cross_Sections::sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  return sigmaCCn(Ebeam, k, p, isProton, 1);
}


double Cross_Sections::sigmaCC2(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  return sigmaCCn(Ebeam, k, p, isProton, 2);
}

double Cross_Sections::sigma_onShell_by_Etheta(double Ebeam, TVector3 k, bool isProton)
{
  double theta=k.Theta();
  double E3 = Ebeam * mN/ (mN + Ebeam*(1.-k.CosTheta()));
  double QSq = 2. * Ebeam * E3 * (1.-k.CosTheta());
  double tau = QSq/(4.*mN*mN);
  double GE = isProton ? GEp(QSq) : GEn(QSq);
  double GM = isProton ? GMp(QSq) : GMn(QSq);
  double epsilon = epsilon = 1./(1.+2.*(1.+tau)*sq(tan(theta/2.)));
  double sigmaMott = cmSqGeVSq * 4* sq(alpha) * E3*E3*E3 * sq(cos(theta/2.))/(Ebeam*QSq*QSq);
  return sigmaMott * (sq(GE) + tau/epsilon * sq(GM))/(1. + tau);
}

double Cross_Sections::GEp(double QSq)
{
  switch (myModel)
    {
    case dipole: 
      return Gdipole(QSq);
    case kelly:
      return Gkelly(QSq,-0.24,10.98,12.82,21.97);
    default:
      std::cerr << "Error in GEp: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double Cross_Sections::GEn(double QSq) // This will use the Galster parameterization
{
  double tau = QSq/(4.*mN*mN);
  return 1.70 * tau / (1. + 3.3 * tau) * Gdipole(QSq); // params from Kelly paper
  //return 1.91 * tau / (1. + 5.6 * tau) * Gdipole(QSq); // I can't find a source for these numbers, so I trust them less.
}

double Cross_Sections::GMp(double QSq)
{
  switch (myModel)
    {
    case dipole: 
      return mu_p * Gdipole(QSq);
    case kelly:
      return mu_p * Gkelly(QSq,0.12,10.97,18.86,6.55);
    default:
      std::cerr << "Error in GMp: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double Cross_Sections::GMn(double QSq)
{
  switch (myModel)
    {
    case dipole: 
      return mu_n * Gdipole(QSq);
    case kelly:
      return mu_n * Gkelly(QSq,2.33,14.72,24.20,84.1);
    default:
      std::cerr << "Error in GMn: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double Cross_Sections::Gdipole(double QSq){ return 1. / sq(1 + QSq/0.71); }

double Cross_Sections::Gkelly(double QSq,double a1, double b1, double b2, double b3)
{
  double tau = QSq/(4.*mN*mN);
  double denom = 1. + b1*tau + b2*tau*tau + b3*tau*tau*tau;
  double numer = 1. + a1*tau;
  return numer/denom;
}

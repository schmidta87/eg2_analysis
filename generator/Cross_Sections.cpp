#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TVector3.h"

#include "Nuclear_Info.h"
#include "Cross_Sections.h"

Cross_Sections::Cross_Sections()
{
}

double Cross_Sections::sq(double x){ return x*x; };
//This is the mostly positvie metric
double Cross_Sections::dot4(double x0, TVector3 x, double y0, TVector3 y)
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
  double GE = (isProton)? Gdipole(QSq) : 1.91 * QSq * Gdipole(QSq) / (4.*sq(mN) + 5.6 * QSq);
  double GM = (isProton)? 2.79*Gdipole(QSq) : -1.91*Gdipole(QSq);
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

  double phi = q.Cross(k).Angle( q.Cross(p) );
  return sigmaMott * ( sq(QSq)/sq(q.Mag2()) * wC +
                       (QSq/(2.*q.Mag2()) + sq(tan(k.Theta()/2.))) * wT +
                       QSq/q.Mag2() * sqrt(QSq/q.Mag2() + sq(tan(k.Theta()/2.))) * wI * cos(phi) +
                       (QSq/q.Mag2() * sq(cos(phi)) + sq(tan(k.Theta()/2.))) * wS
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

double Cross_Sections::Gdipole(double QSq){ return 1. / sq(1 + QSq/0.71); };

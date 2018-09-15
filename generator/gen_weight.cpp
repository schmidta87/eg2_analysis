#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "Nuclear_Info.h"

using namespace std;

// Probability windows
const double Qmin=1.;
const double Qmax=5.;
const double Xmin=1.;
const double Xmax=2.;
const int pCode=2212;
const int nCode=2112;

double sq(double x){ return x*x; };
double sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton);
void do_SXC(int &lead_type, int &rec_type, double r);

int main(int argc, char ** argv)
{
  if (argc !=5)
    {
      cerr << "Wrong number of arguments. Insteady try\n\t"
	   << "gen_weight [A] [Beam energy (GeV)] /path/to/output/file [# of events]\n\n";
      return -1;
    }

  // Read in the arguments
  const double Ebeam=atof(argv[2]);
  const TVector3 v1(0.,0.,Ebeam);
  int nEvents = atoi(argv[4]);
  TFile * outfile = new TFile(argv[3],"RECREATE");

  // Set up the tree
  TTree * outtree = new TTree("genT","Generator Tree");
  Double_t pe[3], q[3], pLead[3], pRec[3], pMiss[3], pCM[3], pRel[3];
  Double_t QSq, xB, nu, pe_Mag, q_Mag, pLead_Mag, pRec_Mag, pMiss_Mag, pCM_Mag, pRel_Mag, theta_pmq, theta_prq, weight;
  Int_t lead_type, rec_type;
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("q",q,"q[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pMiss",pMiss,"pMiss[3]/D");
  outtree->Branch("pCM",pCM,"pCM[3]/D");
  outtree->Branch("pRel",pRel,"pRel[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  outtree->Branch("QSq",&QSq,"QSq/D");
  outtree->Branch("xB",&xB,"xB/D");
  outtree->Branch("nu",&nu,"nu/D");
  outtree->Branch("pe_Mag",&pe_Mag,"pe_Mag/D");
  outtree->Branch("q_Mag",&q_Mag,"q_Mag/D");
  outtree->Branch("pLead_Mag",&pLead_Mag,"pLead_Mag/D");
  outtree->Branch("pRec_Mag",&pRec_Mag,"pRec_Mag/D");
  outtree->Branch("pMiss_Mag",&pMiss_Mag,"pMiss_Mag/D");
  outtree->Branch("pCM_Mag",&pCM_Mag,"pCM_Mag/D");
  outtree->Branch("pRel_Mag",&pRel_Mag,"pRel_Mag/D");
  outtree->Branch("theta_pmq",&theta_pmq,"theta_pmq/D");
  outtree->Branch("theta_prq",&theta_prq,"theta_prq/D");

  // Other chores
  TRandom3 myRand(0);
  Nuclear_Info myInfo(atoi(argv[1]));
  const double mA = myInfo.get_mA();
  const double mAm2 = myInfo.get_mAm2();
  const double sigCM = myInfo.get_sigmaCM();
  
  // Loop over events
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %10000 ==0)
	cerr << "Working on event " << event << "\n";

      // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0
      weight =1.;

      // Decide what kind of proton or neutron pair we are dealing with
      lead_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      rec_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      weight *= 4.;

      // Pick random x, QSq to set up the electron side
      QSq = Qmin + (Qmax-Qmin)*myRand.Rndm();
      xB = Xmin + (Xmax - Xmin)*myRand.Rndm();
      nu = QSq/(2.*mN*xB);
      pe_Mag = Ebeam - nu;
      double cosTheta3 = 1. - QSq/(2.*Ebeam*pe_Mag);
      double phi3 = 2.*M_PI*myRand.Rndm();
      TVector3 v3;
      v3.SetMagThetaPhi(pe_Mag,acos(cosTheta3),phi3);
      pe[0]=v3.X();
      pe[1]=v3.Y();
      pe[2]=v3.Z();
      TVector3 vq = v1 - v3;
      q[0]=vq.X();
      q[1]=vq.Y();
      q[2]=vq.Z();
      q_Mag = vq.Mag();

      // Pick random CM motion
      pCM[0] = myRand.Gaus(0.,sigCM);
      pCM[1] = myRand.Gaus(0.,sigCM);
      pCM[2] = myRand.Gaus(0.,sigCM);
      TVector3 vCM(pCM[0],pCM[1],pCM[2]);
      pCM_Mag = vCM.Mag();
      TVector3 vAm2 = -vCM;
      double EAm2 = sqrt(vCM.Mag2() + sq(mAm2));
      TVector3 vZ = vCM + vq; // 3 momentum of the pair, useful for calculating kinematics
      double X = mA + nu - EAm2;
      double Z = vZ.Mag();
      double YSq = sq(X) - sq(Z);

      // Pick random recoil angles
      double phiRec = myRand.Rndm()*2.*M_PI;
      double cosThetaRec = 2.*(myRand.Rndm()-0.5); // Maybe in the future figure out min and max angles, and restrict
      double thetaRec = acos(cosThetaRec);

      // Determine other angles
      double thetaCM = vCM.Theta();
      double phiCM = vCM.Phi();
      double cosThetaCMRec = sin(thetaCM)*sin(thetaRec) * (cos(phiCM)*cos(phiRec) + sin(phiCM)*sin(phiRec)) + cos(thetaCM)*cosThetaRec;
      double thetaQ = vq.Theta();
      double phiQ = vq.Phi();
      double cosThetaQRec = sin(thetaQ)*sin(thetaRec) * (cos(phiQ)*cos(phiRec) + sin(phiQ)*sin(phiRec)) + cos(thetaQ)*cosThetaRec;
      double thetaZ = vZ.Theta();
      double phiZ = vZ.Phi();
      double cosThetaZRec = sin(thetaZ)*sin(thetaRec) * (cos(phiZ)*cos(phiRec) + sin(phiZ)*sin(phiRec)) + cos(thetaZ)*cosThetaRec;

      // Determine recoil momentum
      double D = sq(X)*(sq(YSq) + sq(2.*mN)*(vZ.Mag2()*sq(cosThetaZRec) - sq(X)));

      if (D<0)
	{
	  weight=0.;
	}
      else
	{
	  double momRec1 = 0.5*(YSq*Z*cosThetaZRec + sqrt(D))/(sq(X) - vZ.Mag2()*sq(cosThetaZRec));
	  double momRec2 = 0.5*(YSq*Z*cosThetaZRec - sqrt(D))/(sq(X) - vZ.Mag2()*sq(cosThetaZRec));
	  if ((momRec1 < 0) && (momRec2 <0))
	    {
	      weight=0.;
	    }
	  else
	    {
	      if (momRec1 < 0)
		pRec_Mag = momRec2;
	      else if (momRec2 < 0)
		pRec_Mag = momRec1;
	      else
		{
		  pRec_Mag = (gRandom->Rndm()>0.5)? momRec1 : momRec2;
		  weight*=2.; // because the solution we picked is half as likely
		}

	      // Define the recoil vector
	      TVector3 vRec;
	      vRec.SetMagThetaPhi(pRec_Mag,acos(cosThetaRec),phiRec);
	      pRec[0]=vRec.X();
	      pRec[1]=vRec.Y();
	      pRec[2]=vRec.Z();

	      // Define the other vectors in the tree
	      pMiss[0] = pCM[0] - pRec[0];
	      pMiss[1] = pCM[1] - pRec[1];
	      pMiss[2] = pCM[2] - pRec[2];
	      pLead[0] = pMiss[0] + q[0];
	      pLead[1] = pMiss[1] + q[1];
	      pLead[2] = pMiss[2] + q[2];
	      pRel[0] = 0.5*(pMiss[0] - pRec[0]);
	      pRel[1] = 0.5*(pMiss[1] - pRec[1]);
	      pRel[2] = 0.5*(pMiss[2] - pRec[2]);
	      TVector3 vLead(pLead[0],pLead[1],pLead[2]);
	      TVector3 vRel(pRel[0],pRel[1],pRel[2]);
	      TVector3 vMiss(pMiss[0],pMiss[1],pMiss[2]);
	      pLead_Mag = vLead.Mag();
	      pMiss_Mag = vMiss.Mag();
	      pRel_Mag = vRel.Mag();
	      theta_pmq = acos((pMiss[0]*q[0] + pMiss[1]*q[1] + pMiss[2]*q[2])/pMiss_Mag /q_Mag);
	      theta_prq = acos((pRec[0]*q[0] + pRec[1]*q[1] + pRec[2]*q[2])/pRec_Mag /q_Mag);

	      double Elead = sqrt(sq(mN) + vLead.Mag2());
	      double Erec = sqrt(sq(mN) + vRec.Mag2());

	      // Calculate the weight
	      weight *= sigmaCC1(Ebeam, v3, vLead, (lead_type==2122)) // eN cross section
		* nu/(2.*xB*Ebeam*pe_Mag) * (Qmax-Qmin) * (Xmax-Xmin) // Jacobian for QSq,xB
		* 1./(4.*sq(M_PI)) // Angular terms
		* ((lead_type==rec_type) ? myInfo.get_pp(pRel_Mag) : myInfo.get_pn(pRel_Mag)) // Contacts
		* vRec.Mag2() * Erec * Elead / fabs(Erec*(pRec_Mag - Z*cosThetaZRec) + Elead*pRec_Mag); // Jacobian for delta fnc.
	    }
	}

      // HERE IS WHERE WE SHOULD DO SINGLE CHARGE EXCHANGE!!!!
      do_SXC(lead_type, rec_type,gRandom->Rndm());

      // Fill the tree
      outtree->Fill();      
    } 	  

  // Clean up
  outtree->Write();
  outfile->Close();
  return 0;
}

double Gdipole(double QSq){ return 1. / sq(1 + QSq/0.71); };

double sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  TVector3 q = TVector3(0.,0.,Ebeam) - k;
  TVector3 pM = p-q;

  double QSq = q.Mag2() - sq(Ebeam - k.Mag());
  double E = sqrt(p.Mag2() + sq(mN));
  double Ebar = sqrt(pM.Mag2() + sq(mN));
  double omegabar = E-Ebar;
  double QSqbar = q.Mag2() - sq(omegabar);

  // Calculate form factors
  double GE = (isProton)? Gdipole(QSq) : 1.91 * QSq * Gdipole(QSq) / (4.*sq(mN) + 5.6 * QSq);
  double GM = (isProton)? 2.79*Gdipole(QSq) : -1.91*Gdipole(QSq);
  double F1 = 0.5 * (GE + QSq*GM/(4.*sq(mN)));
  double kF2 = (GM - GE)/(1. + QSq/(4.*sq(mN)));

  double wC = (sq(E+Ebar)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2)) - q.Mag2()*sq(F1 + kF2))/(4.*E*Ebar);
  double wT = QSqbar*sq(F1 + kF2)/(2.*Ebar*E);
  double wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
  double wI = -p.Mag()*sin(p.Angle(q))*(Ebar + E)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);

  double sigmaMott = cmSqGeVSq * 4. * sq(alpha) * k.Mag2() * sq(cos(k.Theta()/2.)) / sq(QSq);

  double phi = q.Cross(k).Angle( q.Cross(p) );
  return sigmaMott * ( sq(QSq)/q.Mag2() * wC +
                       (QSq/(2.*q.Mag2()) + sq(tan(k.Theta()/2.))) * wT +
                       QSq/q.Mag2() * sqrt(QSq/q.Mag2() + sq(tan(k.Theta()/2.))) * wI * cos(phi) +
                       (QSq/q.Mag2() * sq(cos(phi)) + sq(tan(k.Theta()/2.))) * wS
                       );
}

void do_SXC(int &lead_type, int &rec_type, double r)
{
  const double pPP2PN = 0.048;
  const double pPP2NP = 0.041;
  const double pPP2NN = 0.0029;

  const double pPN2NN = 0.035;
  const double pPN2PP = 0.041;
  const double pPN2NP = 0.0021;

  const double pNP2NN = 0.041;
  const double pNP2PP = 0.035;
  const double pNP2PN = 0.0021;

  const double pNN2PN = 0.041;
  const double pNN2NP = 0.048;
  const double pNN2PP = 0.029;
 
  // Now we do a whole bunch of tests
  if ((lead_type==pCode) && (rec_type==pCode))
    {
      if (r < pPP2NN)
	{
	  lead_type=nCode;
	  rec_type=nCode;
	}
      else if (r < pPP2NN + pPP2NP)
	{
	  lead_type=nCode;
	}
      else if (r < pPP2NN + pPP2NP + pPP2PN)
	{
	  rec_type=nCode;
	}
    }
  else if ((lead_type==nCode) && (rec_type==pCode))
    {
      if (r < pNP2PN)
	{
	  lead_type=pCode;
	  rec_type=nCode;
	}
      else if (r < pNP2PN + pNP2PP)
	{
	  lead_type=pCode;
	}
      else if (r < pNP2PN + pNP2PP + pNP2NN)
	{
	  rec_type=nCode;
	}
    }
  else if ((lead_type==pCode) && (rec_type==nCode))
    {
      if (r < pPN2NP)
	{
	  lead_type=nCode;
	  rec_type=pCode;
	}
      else if (r < pPN2NP + pPN2NN)
	{
	  lead_type=nCode;
	}
      else if (r < pPN2NP + pPN2NN + pPN2PP)
	{
	  rec_type=pCode;
	}
    }
  else if ((lead_type==nCode) && (rec_type==nCode))
    {
      if (r < pNN2PP)
	{
	  lead_type=pCode;
	  rec_type=pCode;
	}
      else if (r < pNN2PP + pNN2PN)
	{
	  lead_type=pCode;
	}
      else if (r < pNN2PP + pNN2PN + pNN2NP)
	{
	  rec_type=pCode;
	}
    }
  else
    {
      cerr << "Invalid nucleon codes. Check and fix. Exiting\n\n\n";
    }
}

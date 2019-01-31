#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <ctype.h>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TVectorT.h"

#include "Nuclear_Info.h"
#include "Cross_Sections.h"
#include "helpers.h"

using namespace std;

// Probability windows
const double Qmin=1.;
const double Qmax=5.;
const double Xmin=1.;
const double Xmax=2.;

const bool doRad=true;

double sigmaCC1(double E1, TVector3 k, TVector3 p, bool isProton);
void do_SXC(int &lead_type, int &rec_type, double r);
double deltaHard(double QSq);

void help_mess()
{
  cerr << "Usage: ./gen_weight [path/to/output.root] [Number of desired events] [optional flags]\n"
       << "Optional flags:\n"
       << "-h: Help\n"
       << "-v: Verbose\n"
       << "-z: Print zero-weight events\n"
       << "-A <Nucleus number>==<12>\n"
       << "-s <Sigma_CM [GeV]>\n"
       << "-C <Nuclear Contact [%]> (Use for Cpp0, Cpn0, Cpn1)\n"
       << "-E <E* [GeV]>\n"
       << "-k <kRel cutoff [GeV]==0.25>\n"
       << "-u <Nuclear potential>==<1> (1=AV18, 2=N2LO, 3=N3LO)\n"
       << "-c <Cross section method>==<cc1>\n"
       << "-f <Form Factor model>==<kelly>\n"
       << "-r: Randomize constants\n"
       << "-p: List output file parameter order\n"
       << "-R: Define contacts only by ratio\n"
       << "-I: Define contacts only by inverted ratio\n";
}

void param_mess()
{
  cerr << "Element | Parameter\n"
       << "---------------------------\n"
       << "      0 | Nucleus\n"
       << "      1 | Sigma_CM [GeV]\n"
       << "      2 | Cpp0 [%]\n"
       << "      3 | Cpn0 [%]\n"
       << "      4 | Cpn1 [%]\n"
       << "      5 | E* [GeV]\n"
       << "      6 | pPP2NP\n"
       << "      7 | pPP2PN\n"
       << "      8 | pPP2NN\n"
       << "      9 | pPN2NN\n"
       << "     10 | pPN2PP\n"
       << "     11 | pPN2NP\n"
       << "     12 | pNP2PP\n"
       << "     13 | pNP2NN\n"
       << "     14 | pNP2PN\n"
       << "     15 | pNN2PN\n"
       << "     16 | pNN2NP\n"
       << "     17 | pNN2PP\n"
       << "     18 | kRel Cut [GeV]\n"
       << "     19 | Nuclear Potential\n"
       << "     20 | Cross Section\n"
       << "     21 | Form Factor Model\n ";
}

int main(int argc, char ** argv)
{
    
  if (argc < 2)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_mess();
      return -1;
    }
  
  if (strcmp(argv[1], "-h")==0)
    {
      help_mess();
      return -1;
    }
  
  if (strcmp(argv[1], "-p")==0)
    {
      param_mess();
      return -1;
    }
  
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Insteady try\n\n";
      help_mess();
      return -1;
    }
  
  // Read in the arguments and flags
  TFile * outfile = new TFile(argv[1],"RECREATE");
  int nEvents = atoi(argv[2]);
  
  bool verbose = false;
  int Anum = 12;
  bool do_sCM = false;
  double sCM;
  bool do_Cs = false;
  std::vector<double> Cs;
  double Estar = 0.;
  double pRel_cut = 0.25;
  int pType = 1;
  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  double rand_flag = false;
  bool byRat = false;
  bool byInv = false;
  bool print_zeros=false;

  int c;  
  while ((c=getopt (argc-2, &argv[2], "hvzA:s:C:E:k:u:f:c:rpRI")) != -1) // First two arguments are not optional flags.
    switch(c)
      {
      case 'h':
	help_mess();
	return -1;
      case 'v':
	verbose = true;
	break;
      case 'z':
	print_zeros = true;
	break;
      case 'A':
	Anum = atoi(optarg);
	break;
      case 's':
	do_sCM = true;
	sCM = atof(optarg);
	break;
      case 'C':
	do_Cs = true;
	Cs.push_back(atof(optarg));
	break;
      case 'E':
	Estar = atof(optarg);
	break;
      case 'k':
        pRel_cut = atof(optarg);
	break;
      case 'u':
	pType = atoi(optarg);
	break;
      case 'c':
	if (strcmp(optarg,"onshell")==0)
	  csMeth=onshell;
	else if (strcmp(optarg,"cc1")==0)
	  csMeth=cc1;
	else if (strcmp(optarg,"cc2")==0)
	  csMeth=cc2;
	else if (atoi(optarg)==1)
	  csMeth=cc1;
	else if (atoi(optarg)==2)
	  csMeth=cc2;
	else
	  {
	    cerr << "Invalid cross section designation. Allowed values are onshell, cc1 and cc2. Aborting...\n";
	    return -1;
	  }
	break;
      case 'f':
	if (strcmp(optarg,"kelly")==0)
	  ffMod=kelly;
	else if (strcmp(optarg,"dipole")==0)
	  ffMod=dipole;
	else{
	  cerr << "Invalid form factor model provided. Allowed options are kelly and dipole. Aborting...\n";
	  return -1;
	}
	break;
      case 'r':
	rand_flag = true;
	break;
      case 'p':
	param_mess();
	return -1;
      case 'R':
	byRat = true;
	break;
      case 'I':
	byInv = true;
	break;
      case '?':
	return -1;
      default:
	abort();
      }


  const double Ebeam=eg2beam;
  const TVector3 v1(0.,0.,Ebeam);
  const double lambda_ei = alpha/M_PI * (log( 4.*Ebeam*Ebeam/(me*me)) - 1.);

  TH1D * h_DeltaEi = new TH1D("deltaEi","ISR;Photon Energy [GeV];Counts",100,0.,0.1);
  TH1D * h_DeltaEf = new TH1D("deltaEf","FSR;Photon Energy [GeV];Counts",100,0.,0.1);

  // Set up the tree
  TTree * outtree = new TTree("genT","Generator Tree");
  Double_t pe[3], q[3], pLead[3], pRec[3], pMiss[3], pCM[3], pRel[3];
  Double_t QSq, xB, nu, pe_Mag, q_Mag, pLead_Mag, pRec_Mag, pMiss_Mag, pCM_Mag, pRel_Mag, theta_pmq, theta_prq, weight;
  Int_t lead_type, rec_type;
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  // Pare down the branches that are not strictly needed
  //outtree->Branch("q",q,"q[3]/D");
  //outtree->Branch("pMiss",pMiss,"pMiss[3]/D");
  //outtree->Branch("pCM",pCM,"pCM[3]/D");
  //outtree->Branch("pRel",pRel,"pRel[3]/D");
  //outtree->Branch("QSq",&QSq,"QSq/D");
  //outtree->Branch("xB",&xB,"xB/D");
  //outtree->Branch("nu",&nu,"nu/D");
  //outtree->Branch("pe_Mag",&pe_Mag,"pe_Mag/D");
  //outtree->Branch("q_Mag",&q_Mag,"q_Mag/D");
  //outtree->Branch("pLead_Mag",&pLead_Mag,"pLead_Mag/D");
  //outtree->Branch("pRec_Mag",&pRec_Mag,"pRec_Mag/D");
  //outtree->Branch("pMiss_Mag",&pMiss_Mag,"pMiss_Mag/D");
  //outtree->Branch("pCM_Mag",&pCM_Mag,"pCM_Mag/D");
  //outtree->Branch("pRel_Mag",&pRel_Mag,"pRel_Mag/D");
  //outtree->Branch("theta_pmq",&theta_pmq,"theta_pmq/D");
  //outtree->Branch("theta_prq",&theta_prq,"theta_prq/D");

  // Other chores
  TRandom3 myRand(0);
  Nuclear_Info myInfo(Anum,pType);
  if (do_sCM)
    myInfo.set_sigmaCM(sCM);
  myInfo.set_Estar(Estar);
  if (byRat)
    myInfo.set_byRatio();
  if (do_Cs)
    { 
      if (byRat)
	{
	  std::vector<double>::size_type isize = 1;
	  std::vector<double>::size_type isize2 = 2;
	  if ((Cs.size() != isize) and (Cs.size() != isize2))
	    {
	      cerr << "If defining by contact ratio, please provide only one value and optionally an uncertainty. Aborting...\n";
	      return 1;
	    }
	  if (Cs.size() == isize)
	    myInfo.set_Ratio(Cs[0]);
	  else
	    myInfo.set_Ratio(Cs[0],Cs[1]);
	}
      else if (byInv)
	{
	  std::vector<double>::size_type isize = 1;
	  std::vector<double>::size_type isize2 = 2;
	  if ((Cs.size() != isize) and (Cs.size() != isize2))
	    {
	      cerr << "If defining by contact ratio, please provide only one value and optionally an uncertainty. Aborting...\n";
	      return 1;
	    }
	  if (Cs.size() == isize)
	    myInfo.set_Ratio_Inverse(Cs[0]);
	  else
	    myInfo.set_Ratio_Inverse(Cs[0],Cs[1]);
	}
      else
	{
	  std::vector<double>::size_type isize = 3;
	  if (Cs.size() != isize)
	    {
	      cerr << "Wrong number of contacts added. Require three values (Cpp0, Cpn0, Cpn1). Aborting...\n";
	      return 1;
	    }
	  myInfo.set_Contacts(Cs[0],Cs[1],Cs[2]);
	}
    }
  if (rand_flag)
    {
    myInfo.randomize();
    pRel_cut = 0.23 + myRand.Uniform()*0.04;
    }

  Cross_Sections myCS(csMeth,ffMod);
  const double mA = myInfo.get_mA();
  const double mAm2 = myInfo.get_mAm2(); // this includes the effect of Estar
  const double sigCM = myInfo.get_sigmaCM();

  // Prepare vector of parameters to be output
  TVectorT<double> params(22);
  params[0] = Anum;
  params[1] = sigCM;
  params[2] = myInfo.get_Cpp0();
  params[3] = myInfo.get_Cpn0();
  params[4] = myInfo.get_Cpn1();
  params[5] = myInfo.get_Estar();
  std::vector<double> Ps = myInfo.get_SCX_Ps();
  params[6] = Ps[0];
  params[7] = Ps[1];
  params[8] = Ps[2];
  params[9] = Ps[3];
  params[10] = Ps[4];
  params[11] = Ps[5];
  params[12] = Ps[6];
  params[13] = Ps[7];
  params[14] = Ps[8];
  params[15] = Ps[9];
  params[16] = Ps[10];
  params[17] = Ps[11];
  params[18] = pRel_cut;
  params[19] = pType;
  params[20] = csMeth;
  params[21] = ffMod;
  
  // Loop over events
  for (int event=0 ; event < nEvents ; event++)
    {
      if ((event %100000 ==0) and verbose)
	cerr << "Working on event " << event << "\n";

      // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0
      weight =1.;

      // Decide what kind of proton or neutron pair we are dealing with
      lead_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      rec_type = (myRand.Rndm() > 0.5) ? pCode:nCode;
      weight *= 4.;

      // Start with the radiation off the incoming electron
      double DeltaEi = doRad ? pow(myRand.Rndm(),1./lambda_ei) * Ebeam : 0.;
      h_DeltaEi->Fill(DeltaEi);
      double Ebeam_eff = Ebeam - DeltaEi;

      // Pick a random xB, QSq prior to FSR
      double xB_eff = Xmin + (Xmax - Xmin)*myRand.Rndm();
      double QSq_eff = Qmin + (Qmax-Qmin)*myRand.Rndm();
      double nu_eff = QSq_eff/(2.*mN*xB_eff);
      double pe_Mag_eff = Ebeam_eff - nu_eff;

       if (pe_Mag_eff < 0.)
	weight=0.;
      else
	{

      // The outgoing electron angle won't change in the peaking approximation
      double cosTheta3 = 1. - QSq_eff/(2.*Ebeam_eff*pe_Mag_eff);
      if (fabs(cosTheta3) > 1.)
	weight=0.;
      else
      {
      
      double phi3 = 2.*M_PI*myRand.Rndm();
      double theta3 = acos(cosTheta3);
      TVector3 v3_eff;
      v3_eff.SetMagThetaPhi(pe_Mag_eff, theta3, phi3);
      TVector3 vq_eff = TVector3(0.,0.,Ebeam_eff) - v3_eff;

      // Sample radiation off the outgoing electron
      double lambda_ef = alpha/M_PI * (log( 4.*pe_Mag_eff*pe_Mag_eff/(me*me)) - 1.);
      double DeltaEf = doRad? pow(myRand.Rndm(),1./lambda_ef) * pe_Mag_eff : 0.;
      h_DeltaEf->Fill(DeltaEf);
      pe_Mag = pe_Mag_eff - DeltaEf;

      // This will allow us to calculate apparent quantities
      QSq = 2. * Ebeam * pe_Mag * (1.-cosTheta3);
      nu = Ebeam - pe_Mag;
      xB = QSq/(2.*mN*nu);

      // Fill into vectors
      TVector3 v3;
      v3.SetMagThetaPhi(pe_Mag,theta3,phi3);
      pe[0]=v3.X();
      pe[1]=v3.Y();
      pe[2]=v3.Z();
      TVector3 vq = v1 - v3;
      q[0]=vq.X();
      q[1]=vq.Y();
      q[2]=vq.Z();
      q_Mag = vq.Mag();

      // Pick random CM motion
      TVector3 vCM_eff(myRand.Gaus(0.,sigCM),myRand.Gaus(0.,sigCM),myRand.Gaus(0.,sigCM));
      TVector3 vAm2 = -vCM_eff;
      double EAm2 = sqrt(vCM_eff.Mag2() + sq(mAm2));
      TVector3 vZ = vCM_eff + vq_eff; // 3 momentum of the pair, useful for calculating kinematics
      double X = mA + nu_eff - EAm2;
      double Z = vZ.Mag();
      double YSq = sq(X) - sq(Z);

      // Pick random recoil angles
      double phiRec = myRand.Rndm()*2.*M_PI;
      double cosThetaRec = 2.*(myRand.Rndm()-0.5); // Maybe in the future figure out min and max angles, and restrict
      double thetaRec = acos(cosThetaRec);

      // Determine other angles
      double thetaCM = vCM_eff.Theta();
      double phiCM = vCM_eff.Phi();
      double cosThetaCMRec = sin(thetaCM)*sin(thetaRec) * (cos(phiCM)*cos(phiRec) + sin(phiCM)*sin(phiRec)) + cos(thetaCM)*cosThetaRec;
      double thetaQ = vq_eff.Theta();
      double phiQ = vq_eff.Phi();
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

	  bool momRec1Valid = (momRec1>=0.) && (X - sqrt(sq(momRec1) + sq(mN)) >=0.) && (sq(X) - sq(Z) + 2.*Z*momRec1*cosThetaZRec > 0.);
	  bool momRec2Valid = (momRec2>=0.) && (X - sqrt(sq(momRec2) + sq(mN)) >=0.) && (sq(X) - sq(Z) + 2.*Z*momRec2*cosThetaZRec >= 0.);

	  if ( (!momRec1Valid) && (!momRec2Valid))
	    {
	      weight=0.;
	    }
	  else
	    {
	      if (!momRec1Valid)
		pRec_Mag = momRec2;
	      else if (!momRec2Valid)
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
	      TVector3 vMiss_eff = vCM_eff - vRec; // This is the true pmiss
	      TVector3 vLead = vMiss_eff + vq_eff; // This is the true plead
	      pLead[0] = vLead.X();
	      pLead[1] = vLead.Y();
	      pLead[2] = vLead.Z();
	      pLead_Mag = vLead.Mag();
	      TVector3 vMiss = vLead - vq; // This is the apparent pmiss
	      pMiss[0] = vMiss.X();
	      pMiss[0] = vMiss.Y();
	      pMiss[0] = vMiss.Z();
	      pMiss_Mag = vMiss.Mag();
	      TVector3 vRel_eff = 0.5*(vMiss_eff - vRec); // This is the true pRel
	      double pRel_eff_Mag = vRel_eff.Mag();

	      // Do a safeguard cut
	      if (pRel_eff_Mag < pRel_cut)
		weight=0.;
	      
	      pRel[0] = 0.5*(pMiss[0] - pRec[0]); // This is the apparent pRel;
	      pRel[1] = 0.5*(pMiss[1] - pRec[1]);
	      pRel[2] = 0.5*(pMiss[2] - pRec[2]);
	      TVector3 vRel(pRel[0],pRel[1],pRel[2]);
	      pRel_Mag = vRel.Mag();

	      pCM[0] = pMiss[0] + pRec[0]; // Apparent pCM
	      pCM[1] = pMiss[1] + pRec[1];
	      pCM[2] = pMiss[2] + pRec[2];
	      TVector3 vCM(pCM[0],pCM[1],pCM[2]);
	      pCM_Mag = vCM.Mag();

	      // These are apparent angles
	      theta_pmq = acos((pMiss[0]*q[0] + pMiss[1]*q[1] + pMiss[2]*q[2])/pMiss_Mag /q_Mag);
	      theta_prq = acos((pRec[0]*q[0] + pRec[1]*q[1] + pRec[2]*q[2])/pRec_Mag /q_Mag);

	      double Elead = sqrt(sq(mN) + vLead.Mag2()); // True values
	      double Erec = sqrt(sq(mN) + vRec.Mag2());

	      // Calculate the weight
	      weight *= myCS.sigma_eN(Ebeam_eff, v3_eff, vLead, (lead_type==pCode)) // eN cross section
		* nu_eff/(2.*xB_eff*Ebeam_eff*pe_Mag_eff) * (Qmax-Qmin) * (Xmax-Xmin) // Jacobian for QSq,xB
		* (doRad ? (1. - deltaHard(QSq_eff)) * pow(Ebeam/sqrt(Ebeam*pe_Mag),lambda_ei) * pow(pe_Mag_eff/sqrt(Ebeam*pe_Mag),lambda_ef) : 1.) // Radiative weights
		* 1./(4.*sq(M_PI)) // Angular terms
		* ((lead_type==rec_type) ? myInfo.get_pp(pRel_eff_Mag) : myInfo.get_pn(pRel_eff_Mag)) // Contacts
		* vRec.Mag2() * Erec * Elead / fabs(Erec*(pRec_Mag - Z*cosThetaZRec) + Elead*pRec_Mag); // Jacobian for delta fnc.
	    }
	}
	}
	}
      
      // HERE IS WHERE WE SHOULD DO SINGLE CHARGE EXCHANGE!!!!
      myInfo.do_SXC(lead_type, rec_type,gRandom->Rndm());
      
      // Fill the tree
      if ((weight > 0.) || print_zeros)
	outtree->Fill();      
    }

  if (verbose)
    {
      cerr << "Listing parameters:\n";
      params.Print();
    }
  
  // Clean up
  params.Write("parameters");
  h_DeltaEi->Write();
  h_DeltaEf->Write();
  outtree->Write();
  outfile->Close();
  return 0;
}

double deltaHard(double QSq)
{
  return 2.*alpha/M_PI * ( -13./12.*log(QSq/(me*me)) + 8./3.);
}

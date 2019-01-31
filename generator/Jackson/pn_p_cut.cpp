#include <cstdlib>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include "AccMap.h"
#include "fiducials.h"
#include "Nuclear_Info.h"
#include "Cross_Sections.h"

using namespace std;

const bool doSmearing=true;
double eSmearing=0.003;
double pSmearing=0.01;
double Tp = 0.53;
double sig_Tp = 0.05;
double Tpp = 0.44;
double sig_Tpp = 0.04;

double sq(double x)
{
  return x*x;
}

bool acceptanceFunction(TVector3 v){

  double theta=v.Theta()*180./TMath::Pi();
  double phi = v.Phi()*180./TMath::Pi();
  double mom=v.Mag();

  double thetaBars[57]={8.6,10.3,11.9,13.6,15.3,17,18.8,20.5,22.3,24,25.8,27.6,29.3,31.1,32.8,34.5,36.2,37.9,39.6,41.2,42.8,44.4,45.9,47.4,49.6,51.9,54.3,56.8,59.4,62,64.7,67.4,70.2,72.9,75.7,78.2,80.8, 
			83.5,86.3,89.3,92.2,95.3,98.4,101.6,104.8,108.8,112,114.9,118.1,121.4,124.7,128,131.4,134.2,136.5,138.8,141}; // cut 5 cm from edge


  double width[57]={32.3,48.1,64,79.8,95.7,106.6,122.4,138.3,154.1,170,185.8,201.7,217.6,233.4,249.3,265.1,281,296.8,312.7,328.5,344.4,360.2,376.1,371.3,378.2,385,391.9,398.7,405.6,412.5,419.3,426.2, 
		    433,439.9,445.1,439.3,433.6,427.8,422,416.3,410.5,404.8,399,393.3,387.5,380.1,363.5,347,330.4,313.9,297.3,280.8,264.2,246.8,235.4,224,212.7};
  double distance[57]={513,509,503,500,498,496,495,494,493,493,493,494,495,496,498,501,504,507,511,515,519,524,529,514,504,495,487,480,473,468,463,460,458,457,457,446,437,428,421,414,409,405,402,400, 
		       399,402,395,389,384,381,378,377,378,379,380,383,385};


  // using cut of Pneutron > 0.35 in the analysis.
  if(mom>-1 && mom<0.35)
    return false;
  if(theta<9 || theta>140)
    return false;

  // put the out of plane angle in the sector region +- 30 deg.
  double PhiOneSector;
  if(TMath::Abs(phi)<30)
    PhiOneSector = phi;
  if(TMath::Abs(phi-60)<30)
    PhiOneSector = phi-60;
  if(TMath::Abs(phi-120)<30)
    PhiOneSector = phi-120;
  if(TMath::Abs(phi-180)<30)
    PhiOneSector = phi - 180;
  if(TMath::Abs(phi+60)<30)
    PhiOneSector = phi+60;
  if(TMath::Abs(phi+120)<30)
    PhiOneSector = phi+120;
  if(TMath::Abs(phi+180)<30)
    PhiOneSector = phi + 180;

  // calculate the minimum and maximam position along each bar.
  double PhiLimitMax[57];
  double PhiLimitMin[57];
  double distFromEdge=10;
  for(int ii=0;ii<57;ii++){
    double tempX = TMath::Sin(thetaBars[ii]*TMath::Pi()/180.) * distance[ii];
    double tempY = width[ii]/2. - distFromEdge;
    PhiLimitMax[ii] = TMath::ATan2(tempY,tempX)*180./TMath::Pi();
    PhiLimitMin[ii] = -1.0 * PhiLimitMax[ii];
  }
  TGraph *g1=new TGraph(57,thetaBars,PhiLimitMax);
  TGraph *g2=new TGraph(57,thetaBars,PhiLimitMin);

  double max = g1->Eval(theta);
  double min = g2->Eval(theta);

  //cout<<" min = "<<min<<" max = "<<max<<" PhiOneSector = "<<PhiOneSector<<endl;
  if(PhiOneSector>min && PhiOneSector<max)
    return true;
  else
    return false;

}

int main(int argc, char ** argv)
{
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tpn_p_cut /path/to/gen/file /path/to/map/file /path/to/out/file\n";
      return -1;
    }

  // Read arguments, create files
  TFile * infile = new TFile(argv[1]);
  AccMap pMap(argv[2]);
  AccMap eMap(argv[2],"e");
  TFile * outfile = new TFile(argv[3],"RECREATE");
  
  // Input Tree
  TTree * inTree = (TTree*)infile->Get("genT");
  Double_t gen_pe[3], gen_pLead[3], gen_pRec[3], gen_weight;
  Int_t lead_type, rec_type;
  inTree->SetBranchAddress("lead_type",&lead_type);
  inTree->SetBranchAddress("rec_type",&rec_type);
  inTree->SetBranchAddress("weight",&gen_weight);
  inTree->SetBranchAddress("pe",gen_pe);
  inTree->SetBranchAddress("pLead",gen_pLead);
  inTree->SetBranchAddress("pRec",gen_pRec);

  // Output Tree
  TTree * T = new TTree("T","Cut Tree");
  Float_t Q2, Xb, Pe[3], Pe_size, theta_e, phi_e, Pp[2][3], Pp_size[2], pq_angle[2], Ep[2], theta_p[2], phi_p[2], nu, q[3];
  Float_t Pmiss_q_angle[2], Pmiss_size[2], Pmiss[2][3];
  Float_t z = 0.;
  Float_t Rp[2][3]={{0.,0.,-22.25},{0.,0.,-22.25}};
  Double_t weightp, weightpn;
  Double_t pRec_Mag;
  Double_t pLead_Mag;
  Int_t nmb = 2;
  int ndet;
  
  T->Branch("Q2",&Q2,"Q2/F");
  T->Branch("Xb",&Xb,"Xb/F");
  T->Branch("nu",&nu,"nu/F");
  T->Branch("q",q,"q[3]/F");
  T->Branch("Pe",Pe,"Pe[3]/F");
  T->Branch("Pe_size",&Pe_size,"Pe_size/F");
  T->Branch("theta_e",&theta_e,"theta_e/F");
  T->Branch("phi_e",&phi_e,"phi_e/F");
  T->Branch("Pp",Pp,"Pp[2][3]/F");
  T->Branch("Pp_size",Pp_size,"Pp_size[2]/F");
  T->Branch("pq_angle",pq_angle,"pq_angle[2]/F");
  T->Branch("Ep",Ep,"Ep[2]/F");
  T->Branch("theta_p",theta_p,"theta_p[2]/F");
  T->Branch("phi_p",phi_p,"phi_p[2]/F");
  T->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[2]/F");
  T->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[2]/F");
  T->Branch("Pmiss",Pmiss,"Pmiss[2][3]/F");
  T->Branch("Rp",Rp,"Rp[2][3]/F");
  T->Branch("weightp",&weightp,"weightp/D");
  T->Branch("weightpn",&weightpn,"weightpn/D");
  T->Branch("pRec",&pRec_Mag,"pRec/D");
  T->Branch("pLead",&pLead_Mag,"pLead/D");

  TRandom3 myRand(0);
  double a = 0.03;
  double b = 0.06;
  Cross_Sections myCS; 
  
  // Loop over all events
  const int nEvents = inTree->GetEntries(); // this is a key number for the weight
  for (int event=0 ; event < nEvents ; event++)
    {
      if (event %100000==0)
	cerr << "Working on event " << event << " out of " << nEvents <<"\n";
      
      inTree->GetEvent(event);
      
      // Require a lead proton
      if (lead_type != pCode)
	continue;

      // Create vectors for the particles
      TVector3 vrec(gen_pRec[0],gen_pRec[1],gen_pRec[2]);
      TVector3 ve(gen_pe[0],gen_pe[1],gen_pe[2]);
      TVector3 vq = TVector3(0.,0.,eg2beam) - ve;
      TVector3 vlead(gen_pLead[0],gen_pLead[1],gen_pLead[2]);
      
      // Smearing

      if (doSmearing)
	{
	  ve *= (1. + eSmearing * myRand.Gaus() * ve.Mag());
	  gen_pe[0] = ve.X();
	  gen_pe[1] = ve.Y();
	  gen_pe[2] = ve.Z();

	  vlead *= (1. + pSmearing * myRand.Gaus() * vlead.Mag());
	  gen_pLead[0] = vlead.X();
	  gen_pLead[1] = vlead.Y();
	  gen_pLead[2] = vlead.Z();

	  vrec *= (1. + pSmearing * myRand.Gaus() * vrec.Mag());
	  gen_pRec[0] = vrec.X();
	  gen_pRec[1] = vrec.Y();
	  gen_pRec[2] = vrec.Z();
	}

      // Fiducial cuts
      if (!accept_electron(ve))
	continue;
      if (!accept_proton(vlead))
	continue;

      // Require a recoil neutron, with fiducial cuts
      ndet = 1;
      
      if (rec_type == pCode)
	ndet = 0;
      else if (!acceptanceFunction(vrec))
	ndet = 0;

      if (gen_weight <= 0.)
	continue;

      // Load up tree
      TVector3 vmiss=vlead-vq;
      TVector3 vcm=vmiss+vrec;
      TVector3 vrel=0.5*(vmiss-vrec);
      
      double gen_pMiss_Mag = vmiss.Mag();
      double gen_pe_Mag = ve.Mag();
      double gen_nu = eg2beam - ve.Mag();
      double gen_QSq = 2. * eg2beam * gen_pe_Mag * (1. - ve.CosTheta());
      double gen_xB = gen_QSq/(2.*mN*gen_nu);
      double gen_q_Mag = vq.Mag();
      double gen_pLead_Mag = vlead.Mag();
      double gen_pRec_Mag = vrec.Mag();
      double gen_ELead = sqrt(gen_pLead_Mag*gen_pLead_Mag+mN*mN);
      double m_miss = sqrt((gen_nu+2*mN-gen_ELead)*(gen_nu+2*mN-gen_ELead)-gen_pMiss_Mag*gen_pMiss_Mag);
      
      // Kinematical cuts
      if (gen_xB < 1.2)
	continue;
      if (m_miss > 1.1)
	continue;
      if (0.62 > gen_pLead_Mag/gen_q_Mag)
	continue;
      if (gen_pLead_Mag/gen_q_Mag > 0.96)
	continue;
      if (vlead.Angle(vq)  > 25.*M_PI/180.)
	continue;

      // Apply weight for cross sections
      double gen_theta = ve.Theta();
      double tau = gen_QSq/(4*mN*mN);
      double epsilon = 1/(1.0+2.0*(1.+tau)*sq(tan(gen_theta/2)));

      weightp = gen_weight * pMap.accept(vlead) * Tp * 1.E33;
      weightpn = gen_weight * pMap.accept(vlead) * Tpp * 1.E33 * ndet;

      // Load up tree
      pRec_Mag = gen_pRec_Mag;
      pLead_Mag = gen_pLead_Mag;

      // Load up tree
      Q2 = gen_QSq;
      Xb = gen_xB;
      nu = gen_nu;
      Pe_size = gen_pe_Mag;
      theta_e = ve.Theta() * 180./M_PI;
      phi_e = ve.Phi()*180./M_PI;
      if (phi_e < -30.) phi_e += 360;
      Pp_size[0] = gen_pLead_Mag;
      Pp_size[1] = gen_pRec_Mag;
      pq_angle[0] = vq.Angle(vlead)*180./M_PI;
      pq_angle[1] = vq.Angle(vrec)*180./M_PI;
      Ep[0] = sqrt(gen_pLead_Mag*gen_pLead_Mag + mN*mN);
      Ep[1] = sqrt(gen_pRec_Mag*gen_pRec_Mag + mN*mN);
      theta_p[0] = vlead.Theta()*180./M_PI;
      theta_p[0] = vrec.Theta()*180./M_PI;
      phi_p[0] = vlead.Phi()*180./M_PI;
      if (phi_p[0] < -30.) phi_p[0] += 360;
      phi_p[1] = vrec.Phi()*180./M_PI;
      if (phi_p[1] < -30.) phi_p[1] += 360;
      for (int i=0 ; i<3 ; i++)
	{
	  q[i] = vq[i];
	  Pe[i] = gen_pe[i];
	  Pp[0][i] = vlead[i];
	  Pp[1][i] = vrec[i];
	  Pmiss[0][i] = vlead[i] - vq[i];
	  Pmiss[1][i] = vrec[i] - vq[i];
	}

      Pmiss_q_angle[0] = (vlead - vq).Angle(vq) * 180./M_PI;
      Pmiss_q_angle[1] = (vrec - vq).Angle(vq) * 180./M_PI;
      
      Pmiss_size[0] = (vlead-vq).Mag();
      Pmiss_size[1] = (vrec - vq).Mag();
      
      T->Fill();
    }
  
  infile->Close();
  
  outfile->cd();
  T->Write();
  outfile->Close();

  return 0;
}


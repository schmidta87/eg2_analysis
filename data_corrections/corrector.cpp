#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include "constants.h"

using namespace std;

double correct_theta(TVector3 v_uncorr);

// Subtract these from the measured electron momentum for each sector.
const double Ee_offsets[6]={-0.00597309,0.0211157,-0.00535153,-0.0260952,-0.00339546,0.0095037};

int main(int argc, char ** argv)
{
  // Must take in an uncorrected file, and spit out a corrected file.
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "\tcorrector /path/to/uncorrected/data/file /path/to/corrected/output/file\n\n";
      return -1;
    }

  // Get the in tree
  TFile * inFile = new TFile(argv[1]);
  TTree * inTree = (TTree*)inFile->Get("T");
  const int maxNP = 2;

  // Memory to store the data
  Int_t nmb, target_type, ID[maxNP];
  Float_t Q2, Xb, Pe[3], Pe_size, Ee, phi_e, theta_e, EC_in, EC_out, EC_tot,
    electron_z, Pp[maxNP][3],Pp_size[maxNP],pq_angle[maxNP],Ep[maxNP], Rp[maxNP][3],
    phi_p[maxNP], theta_p[maxNP], TimeCorr4[maxNP], SC_Edep[maxNP], Nu, q[3],
    q_size, theta_q, phi_q, theta_Pmiss[maxNP], phi_Pmiss[maxNP], Pmiss[maxNP][3],
    Pmiss_size[maxNP], Pmiss_q_angle[maxNP];

  // Set up the output tree
  TFile * outFile = new TFile(argv[2],"RECREATE");
  TTree * outTree = new TTree("T",inTree->GetTitle());
  
  // First, connect the input and output trees for data items with no correction
  inTree->SetBranchAddress("nmb",&nmb);
  outTree->Branch("nmb",&nmb,"nmb/I");
  inTree->SetBranchAddress("target_type",&target_type);
  outTree->Branch("target_type",&target_type,"target_type/I");
  inTree->SetBranchAddress("ID",ID);
  outTree->Branch("ID",&ID,"ID[nmb]/I");
  inTree->SetBranchAddress("EC_in",&EC_in);
  outTree->Branch("EC_in",&EC_in,"EC_in/F");
  inTree->SetBranchAddress("EC_out",&EC_out);
  outTree->Branch("EC_out",&EC_out,"EC_out/F");
  inTree->SetBranchAddress("EC_tot",&EC_tot);
  outTree->Branch("EC_tot",&EC_tot,"EC_tot/F");
  inTree->SetBranchAddress("electron_z",&electron_z);
  outTree->Branch("electron_z",&electron_z,"electron_z/F");
  inTree->SetBranchAddress("Rp",Rp);
  outTree->Branch("Rp",Rp,"Rp[nmb][3]/F");
  inTree->SetBranchAddress("TimeCorr4",TimeCorr4);
  outTree->Branch("TimeCorr4",TimeCorr4,"TimeCorr4[nmb]/F");
  inTree->SetBranchAddress("SC_Edep",SC_Edep);
  outTree->Branch("SC_Edep",SC_Edep,"SC_Edep[nmb]/F");

  // Next, hook up the output tree for everything else
  outTree->Branch("Q2",&Q2,"Q2/F");
  outTree->Branch("Xb",&Xb,"Xb/F");
  outTree->Branch("Pe_size",&Pe_size,"Pe_size/F");
  outTree->Branch("Ee",&Ee,"Ee/F");
  outTree->Branch("phi_e",&phi_e,"phi_e/F");
  outTree->Branch("theta_e",&theta_e,"theta_e/F");
  outTree->Branch("Nu",&Nu,"Nu/F");
  outTree->Branch("q_size",&q_size,"q_size/F");
  outTree->Branch("theta_q",&theta_q,"theta_q/F");
  outTree->Branch("phi_q",&phi_q,"phi_q/F");
  outTree->Branch("Pe",Pe,"Pe[3]/F");
  outTree->Branch("q",q,"q[3]/F");
  outTree->Branch("Pp",Pp,"Pp[nmb][3]/F");
  outTree->Branch("Pp_size",Pp_size,"Pp_size[nmb]/F");
  outTree->Branch("pq_angle",pq_angle,"pq_angle[nmb]/F");
  outTree->Branch("Ep",Ep,"Ep[nmb]/F");
  outTree->Branch("phi_p",phi_p,"phi_p[nmb]/F");
  outTree->Branch("theta_p",theta_p,"theta_p[nmb]/F");
  outTree->Branch("theta_Pmiss",theta_Pmiss,"theta_Pmiss[nmb]/F");
  outTree->Branch("phi_Pmiss",phi_Pmiss,"phi_Pmiss[nmb]/F");
  outTree->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[nmb]/F");
  outTree->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[nmb]/F");
  outTree->Branch("Pmiss",Pmiss,"Pmiss[nmb][3]/F");

  // Lastly, hook up the minimal set of correctable variables
  inTree->SetBranchAddress("Pe",Pe);
  inTree->SetBranchAddress("Pp",Pp);

  // Loop over events
  TVector3 vps_uncorr[maxNP];
  TVector3 vps[maxNP];
  for (int event=0 ; event < inTree->GetEntries() ; event++)
    {
      // Load the input vectors
      inTree->GetEvent(event);

      // Create uncorrected vectors
      TVector3 ve_uncorr(Pe[0],Pe[1],Pe[2]);
      for (int p=0 ; p< nmb ; p++)
	vps_uncorr[p].SetXYZ(Pp[p][0],Pp[p][1],Pp[p][2]);

      // Determine the sector
      double temp_phi=ve_uncorr.Phi();
      if (temp_phi < -M_PI/6.) temp_phi+= 2.*M_PI;
      int sector = temp_phi/(M_PI/3.);

      // Do the correction
      TVector3 ve;
      ve.SetMagThetaPhi(ve_uncorr.Mag()-Ee_offsets[sector],correct_theta(ve_uncorr), ve_uncorr.Phi());
      
      for (int p=0 ; p<nmb ; p++)
	  vps[p].SetMagThetaPhi(vps_uncorr[p].Mag(),correct_theta(vps_uncorr[p]),vps_uncorr[p].Phi());

      // Produced derived variables
      TVector3 vq = TVector3(0.,0.,eg2beam) - ve;
     
      // Load up all of the output tree variables
      Pe_size = ve.Mag();
      Ee = ve.Mag();
      theta_e = ve.Theta() * 180./M_PI;
      Nu = eg2beam - Ee;
      q_size = vq.Mag();
      theta_q = vq.Theta() * 180./M_PI;
      Q2 = vq.Mag2() - Nu*Nu;
      Xb = Q2 / (2. * mN * Nu);

      phi_e = ve.Phi() * 180./M_PI;
      if (phi_e < -30) phi_e += 360.;

      phi_q = vq.Phi() * 180./M_PI;
      if (phi_q < -30) phi_q += 360.;

      // Those that are length 3
      for (int i=0 ; i<3 ; i++)
	{
	  Pe[i] = ve[i];
	  q[i] = vq[i];
	}

      // Those that depend on specific proton
      for (int p=0 ; p<nmb ; p++)
	{
	  Pp_size[p] = vps[p].Mag();
	  pq_angle[p] = vps[p].Angle(vq) * 180./M_PI;
	  Ep[p] = sqrt( vps[p].Mag2() + mN*mN);
	  theta_p[p] = vps[p].Theta() * 180./M_PI;

	  phi_p[p] = vps[p].Phi() * 180./M_PI;
	  if (phi_p[p] < -30.) phi_p[p] += 360.;

	  TVector3 vmiss = vps[p] - vq;

	  Pmiss_size[p] = vmiss.Mag();
	  theta_Pmiss[p] = vmiss.Theta() * 180./M_PI;
	  Pmiss_q_angle[p] = vmiss.Angle(q) * 180./M_PI;
	  
	  phi_Pmiss[p] = vmiss.Phi() * 180./M_PI;
	  if (phi_Pmiss[p] < -30.) phi_Pmiss[p] += 360.;

	  for (int i=0 ; i<3 ; i++)
	    {
	      Pp[p][i] = vps[p][i];
	      Pmiss[p][i] = vmiss[i];
	    }
	}

      outTree->Fill();
    }

  inFile->Close();
  outFile->cd();
  outTree->Write();
  outFile->Close();
  return 0;
}

double correct_theta(TVector3 v_uncorr)
{  
  // Key parameters (all angles need to be in degrees)
  const double ThetaCorrection[6][4] = {{-6.898465e-03,-5.093691e-04,-4.437409e-05,-1.085586e-06},
					{-3.376092e-01,1.085781e-02,-9.263794e-05,4.445614e-08},
					{-6.039766e-01,6.016822e-03,4.450338e-05,-4.373013e-07},
					{-3.417601e-01,2.250364e-03,1.128394e-05,-7.500616e-08},
					{-1.426398e+00,5.806233e-03,2.420251e-05,-9.972778e-08},
					{-3.954411e+00,1.404261e-02,4.487579e-05,-1.591096e-07}};

  double phi = v_uncorr.Phi();
  if (phi < -M_PI/6.) phi += 2.*M_PI;
  double phi_deg = phi*180./M_PI;
  
  int n_sector = phi_deg/60.;
  
  double theta = v_uncorr.Theta();

  double delta_theta = ( ThetaCorrection[n_sector][0] +
			 ThetaCorrection[n_sector][1]*pow(phi_deg,1) +
			 ThetaCorrection[n_sector][2]*pow(phi_deg,2) +
			 ThetaCorrection[n_sector][3]*pow(phi_deg,3)) * M_PI/180.;
  
  return theta + delta_theta;
}


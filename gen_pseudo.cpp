#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorT.h"
#include "TGraphAsymmErrors.h"

#include "AccMap.h"
#include "constants.h"

using namespace std;

// Global variables for the parameters
int n_trials = 200.;

int main(int argc, char ** argv)
{
  if (argc != 11)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "/path/to/1p/file /path/to/2p/file /path/to/map/file /path/to/output/file [# iter] [a1] [a2] [b1] [b2] [sigPerp]\n\n";
      return -1;
    }

  // Get the important info out of the arguments
  const int samples = atoi(argv[5]);
  const double a1 = atof(argv[6]);
  const double a2 = atof(argv[7]);
  const double b1 = atof(argv[8]);
  const double b2 = atof(argv[9]);
  const double sigPerp = atof(argv[10]);
  TFile * outfile = new TFile(argv[4],"RECREATE");
  AccMap myMap(argv[3]);

  // Vectors which describe all of the useful info for the e'p events
  vector<TVector3> ep_pmiss_list;
  vector<TVector3> ep_q_list;
  vector<TVector3> ep_lon_list;
  vector<TVector3> ep_inp_list;
  vector<TVector3> ep_oop_list;
  
  // Info for the e'pp events
  double * epp_acc_weights;
  double * epp_pcm_lon_list;
  double * epp_pcm_inp_list;
  double * epp_pcm_oop_list;

  // Initialize some histograms
  TH2D * hist_ep_events = new TH2D("ep","ep events;p_miss [GeV];angle between p_miss and q [deg];Counts",n_pmiss_bins,0.3,1.,30,90.,180.);
  TH2D * hist_epp_cm_lon = new TH2D("cm_lon","epp events;p_miss [GeV];p_cm_lon [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_inp = new TH2D("cm_inp","epp events;p_miss [GeV];p_cm_inp [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_oop = new TH2D("cm_oop","epp events;p_miss [GeV];p_cm_oop [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_model_cm_lon = new TH2D("model_cm_lon","Underlying model likelihood;p_miss [GeV];p_cm_lon [GeV];Counts",7,0.3,1.,120,-1.2,1.2);
  TH2D * hist_model_cm_inp = new TH2D("model_cm_inp","Underlying model likelihood;p_miss [GeV];p_cm_inp [GeV];Counts",7,0.3,1.,120,-1.2,1.2);
  TH2D * hist_model_cm_oop = new TH2D("model_cm_oop","Underlying model likelihood;p_miss [GeV];p_cm_oop [GeV];Counts",7,0.3,1.,120,-1.2,1.2);
  TH2D * hist_model_fine_cm_lon = new TH2D("model_fine_cm_lon","Underlying model likelihood;p_miss [GeV];p_cm_lon [GeV];Counts",n_pmiss_bins,0.3,1.,120,-1.2,1.2);
  TH2D * hist_model_fine_cm_inp = new TH2D("model_fine_cm_inp","Underlying model likelihood;p_miss [GeV];p_cm_inp [GeV];Counts",n_pmiss_bins,0.3,1.,120,-1.2,1.2);
  TH2D * hist_model_fine_cm_oop = new TH2D("model_fine_cm_oop","Underlying model likelihood;p_miss [GeV];p_cm_oop [GeV];Counts",n_pmiss_bins,0.3,1.,120,-1.2,1.2);
  TH1D * hist_pmiss_gen = new TH1D("pseudo_pmiss_gen","Generated pseudo events;p_miss [GeV];Counts",n_pmiss_bins,0.3,1.);
  TH1D * hist_pmiss_acc = new TH1D("pseudo_pmiss_acc","Accepted pseudo events;p_miss [GeV];Counts",n_pmiss_bins,0.3,1.);
  hist_pmiss_gen->Sumw2();
  hist_pmiss_acc->Sumw2();
  TGraphAsymmErrors * pp2p_corr = new TGraphAsymmErrors();
  pp2p_corr->SetName("pseudo_pmiss_corr");
  pp2p_corr->SetTitle("pp/p correction for pseudo data;p_miss [GeV];Correction");

  // Initialize random guess
  TRandom3 myRand(0);

  // Read in e'p events
  cerr << "Loading 1p data file " << argv[1] << " ...\n";

  TFile * fep = new TFile(argv[1]);
  TTree * tep = (TTree*)fep->Get("T");
  Float_t pPVec[2][3]; // Proton momenta
  Float_t qVec[3]; // q
  Float_t pMissVec[2][3];
  Float_t vertices[2][3]; // position 3 vectors of proton vertices
  tep->SetBranchAddress("q",qVec);
  tep->SetBranchAddress("Pmiss",pMissVec);
  tep->SetBranchAddress("Rp",vertices);
  tep->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tep->GetEntries() ; i++)
    {
      tep->GetEvent(i);

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);
      
      // Make Or's cuts on vertex position
      if(!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
	continue;

      ep_pmiss_list.push_back(pmiss);
      ep_q_list.push_back(q);

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = oop.Cross(lon).Unit();
      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);

      hist_ep_events->Fill(pmiss.Mag(),pmiss.Angle(q)*180./M_PI);
    }
  fep->Close();

  // Read in epp file
  cerr << "Loading 2p data file " << argv[2] << " ...\n";

  TFile *fepp = new TFile(argv[2]);
  TTree *tepp = (TTree*)fepp->Get("T");
  tepp->SetBranchAddress("q",qVec);
  tepp->SetBranchAddress("Pmiss",pMissVec);
  tepp->SetBranchAddress("Rp",vertices);
  tepp->SetBranchAddress("Pp",pPVec);
  for (int i=0 ; i<tepp->GetEntries() ; i++)
    {
      tepp->GetEvent(i);  

      TVector3 plead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);

      // Make Or's cut on vertex position
      if (!((fabs(vertices[0][2]+22.25)<2.25) && (plead.Mag() < 2.4)))
	continue;

      ep_pmiss_list.push_back(pmiss);
      ep_q_list.push_back(q);
      double pmiss_mag = pmiss.Mag();

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = oop.Cross(lon).Unit();

      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);

      hist_ep_events->Fill(pmiss_mag,pmiss.Angle(q)*180./M_PI);
    }
  fepp->Close();

  int n_ep_events = ep_pmiss_list.size();

  cerr << "Done reading input data.\n";
  cerr << "Read in " << n_ep_events << " ep events.\n";

  // Create memory for the acceptance weights
  epp_acc_weights = new double[n_trials * n_ep_events ];
  epp_pcm_lon_list = new double[n_trials * n_ep_events ];
  epp_pcm_inp_list = new double[n_trials * n_ep_events ];
  epp_pcm_oop_list = new double[n_trials * n_ep_events ];
  double weight_sum=0.;

  // Loop over the e'p events and create pseudodata
  for (int i=0 ; i < n_ep_events; i++)
    {
      double pmiss = ep_pmiss_list[i].Mag();
      double sigLong = a1 * (pmiss - 0.6) + a2;
      double muLong = b1 * (pmiss - 0.6) + b2;

      for (int j=0 ; j<n_trials ; j++)
	{
	  double temp_pcm_lon = sigLong * myRand.Gaus() + muLong;
          double temp_pcm_inp = sigPerp * myRand.Gaus();
          double temp_pcm_oop = sigPerp * myRand.Gaus();

	  TVector3 temp_pcm = temp_pcm_lon * ep_lon_list[i] + temp_pcm_inp * ep_inp_list[i] + temp_pcm_oop * ep_oop_list[i];
          TVector3 temp_prec = temp_pcm - ep_pmiss_list[i];

	  double weight = myMap.recoil_accept(temp_prec);
	  hist_pmiss_acc->Fill(pmiss,weight);
	  hist_pmiss_gen->Fill(pmiss,1.);
	    
	  //cout << pmiss << " " << temp_prec.Mag() << " " << temp_pcm_lon << " " << temp_pcm_inp << " " << temp_pcm_oop << " " << weight << "\n";

	  weight_sum += weight;
          epp_acc_weights[i*n_trials + j] = weight_sum;
          epp_pcm_lon_list[i*n_trials + j] = temp_pcm_lon;
          epp_pcm_inp_list[i*n_trials + j] = temp_pcm_inp;
          epp_pcm_oop_list[i*n_trials + j] = temp_pcm_oop;	  

	  hist_model_cm_lon->Fill(pmiss,temp_pcm_lon,weight);
	  hist_model_cm_inp->Fill(pmiss,temp_pcm_inp,weight);
	  hist_model_cm_oop->Fill(pmiss,temp_pcm_oop,weight);

	  hist_model_fine_cm_lon->Fill(pmiss,temp_pcm_lon,weight);
	  hist_model_fine_cm_inp->Fill(pmiss,temp_pcm_inp,weight);
	  hist_model_fine_cm_oop->Fill(pmiss,temp_pcm_oop,weight);
	}
    }

  // Renormalize the model histograms
  for (int i=1 ; i<=n_pmiss_bins ; i++)
    {
      double temp_integral=hist_pmiss_acc->GetBinContent(i);
      for (int j=1 ; j<120 ; j++)
	{
	  int bin = hist_model_fine_cm_lon->GetBin(i,j);
	  hist_model_fine_cm_lon->SetBinContent(bin,hist_model_fine_cm_lon->GetBinContent(bin)/temp_integral);
	  hist_model_fine_cm_inp->SetBinContent(bin,hist_model_fine_cm_inp->GetBinContent(bin)/temp_integral);
	  hist_model_fine_cm_oop->SetBinContent(bin,hist_model_fine_cm_oop->GetBinContent(bin)/temp_integral);
	}
    }

  // Calculate the model's pp to p acceptance correction
  pp2p_corr->BayesDivide(hist_pmiss_acc,hist_pmiss_gen);

  // Renormalize acceptance weights
  for (int i=0 ; i < n_ep_events*n_trials ; i++)
    epp_acc_weights[i] *= 1./weight_sum;

  // Prep output tree
  TTree * outtree = new TTree("T","Pseudodata");
  outtree->Branch("Pp",pPVec,"Pp[2][3]/F");
  outtree->Branch("Rp",vertices,"Rp[2][3]/F");
  outtree->Branch("Pmiss",pMissVec,"Pmiss[2][3]/F");
  outtree->Branch("q",qVec,"q[3]/F");

  // Pick random epp_events  
  for (int i=0 ; i < samples ; i++)
    {
      double r = myRand.Rndm();
      
      int index;
      for (index =0 ; index < n_ep_events*n_trials ; index++)
        if (r < epp_acc_weights[index])
          break;
      
      int ep_event = index / n_trials;

      TVector3 pLead = ep_pmiss_list[ep_event] + ep_q_list[ep_event];
      TVector3 pCM = epp_pcm_lon_list[index]*ep_lon_list[ep_event] + 
	epp_pcm_inp_list[index]*ep_inp_list[ep_event] + 
	epp_pcm_oop_list[index]*ep_oop_list[ep_event];
      TVector3 pRec = pCM - ep_pmiss_list[ep_event];
	
      // Load the important data
      qVec[0] = ep_q_list[ep_event].X();
      qVec[1] = ep_q_list[ep_event].Y();
      qVec[2] = ep_q_list[ep_event].Z();
      vertices[0][0] = vertices[0][1] = vertices[0][2] = vertices[1][0] = vertices[1][1] = vertices[1][2] = -22.5;
      pPVec[0][0] = pLead.X();
      pPVec[0][1] = pLead.Y();
      pPVec[0][2] = pLead.Z();
      pPVec[1][0] = pRec.X();
      pPVec[1][1] = pRec.Y();
      pPVec[1][2] = pRec.Z();
      pMissVec[0][0] = ep_pmiss_list[ep_event].X();
      pMissVec[0][1] = ep_pmiss_list[ep_event].Y();
      pMissVec[0][2] = ep_pmiss_list[ep_event].Z();

      // Fill some histograms
      hist_epp_cm_lon->Fill(ep_pmiss_list[ep_event].Mag(),epp_pcm_lon_list[index]);
      hist_epp_cm_inp->Fill(ep_pmiss_list[ep_event].Mag(),epp_pcm_inp_list[index]);
      hist_epp_cm_oop->Fill(ep_pmiss_list[ep_event].Mag(),epp_pcm_oop_list[index]);
      
      outtree->Fill();
    }

  // Write out histograms
  outfile->cd();
  hist_ep_events->Write();
  hist_epp_cm_lon->Write();
  hist_epp_cm_inp->Write();
  hist_epp_cm_oop->Write();
  hist_model_cm_lon->Write();
  hist_model_cm_inp->Write();
  hist_model_cm_oop->Write();
  hist_model_fine_cm_lon->Write();
  hist_model_fine_cm_inp->Write();
  hist_model_fine_cm_oop->Write();
  hist_pmiss_gen->Write();
  hist_pmiss_acc->Write();
  pp2p_corr->Write();
  outtree->Write();
  outfile->Close();

  // Clean up
  delete[] epp_acc_weights;
  delete[] epp_pcm_lon_list;
  delete[] epp_pcm_inp_list;
  delete[] epp_pcm_oop_list;

  return 0;
}

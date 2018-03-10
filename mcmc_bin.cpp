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

#include "constants.h"
#include "fiducials.h"
#include "AccMap.h"

using namespace std;

// Useful math
double logGaus(double x, double s, double m);

// Global variables for the parameters
int n_epp_events; 
int n_ep_events;
int n_ep_events_per_slice[n_slices];
int n_pseudo_per_ep[n_slices];
int this_sample;
double * current_params;
double * new_params;
double * random_gauss;
long double bestPosterior;
TH2D * bestLonModel;
TH2D * bestInpModel;
TH2D * bestOopModel;

// Vectors which describe all of the useful info for the e'p events
vector<TVector3> ep_pmiss_list;
vector<TVector3> ep_lon_list;
vector<TVector3> ep_inp_list;
vector<TVector3> ep_oop_list;

// Vectors which describe all of the useful info for the e'pp events
vector <double> epp_pmiss_list;
vector <double> epp_pcm_lon_list;
vector <double> epp_pcm_oop_list;
vector <double> epp_pcm_inp_list;

// Acceptance Map
AccMap * myMap;

// Random number generator
TRandom3 * prand;

double jumpscale;

//Box size for the five-dimensional space
const double min_a1 = -0.2;
const double min_a2 = 0.;
const double min_b1 = -0.5;
const double min_b2 = 0.;
const double min_sigPerp = 0.;
const double max_a1 = 0.8;
const double max_a2 = 0.4;
const double max_b1 = 1.5;
const double max_b2 = 0.5;
const double max_sigPerp = 0.4;
const double pcm_bound = 1.;

// Helper functions
long double posterior(double * param_set);
void pick_new_step();
void initialize_params();

int main(int argc, char ** argv)
{
  if (argc != 7)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "/path/to/1p/file /path/to/2p/file /path/to/map/file /path/to/output/file [# iter] [jump-scale]\n\n";
      return -1;
    }

  // Read in the important info
  TFile * outfile = new TFile(argv[4],"RECREATE");
  int samples = atoi(argv[5]);
  jumpscale = atoi(argv[6]);

  // Initialize some histograms
  TH2D * hist_ep_events = new TH2D("ep","ep events;p_miss [GeV];angle between p_miss and q [deg];Counts",n_slices,0.3,1.,30,90.,180.);
  TH2D * hist_epp_events = new TH2D("epp","epp events;p_miss [GeV];angle between p_miss and q [deg];Counts",n_slices,0.3,1.,30,90.,180.);
  TH2D * hist_epp_recoils = new TH2D("epp_rec","epp events;p_rec [GeV];angle between p_rec and q [deg];Counts",n_slices,0.3,1.,30,0.,180.);
  TH2D * hist_epp_cm_lon = new TH2D("cm_lon","epp events;p_miss [GeV];p_cm_lon [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_inp = new TH2D("cm_inp","epp events;p_miss [GeV];p_cm_inp [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  TH2D * hist_epp_cm_oop = new TH2D("cm_oop","epp events;p_miss [GeV];p_cm_oop [GeV];Counts",7,0.3,1.,24,-1.2,1.2);
  bestPosterior = 0.;
  bestLonModel = NULL;
  bestInpModel = NULL;
  bestOopModel = NULL;

  // Initialize random guess
  prand = new TRandom3(0);

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
      
      TVector3 pLead(pPVec[0][0],pPVec[0][1], pPVec[0][2]);      
      TVector3 pmiss(pMissVec[0][0],pMissVec[0][1],pMissVec[0][2]);
      TVector3 q(qVec[0],qVec[1],qVec[2]);

      // Make Or's cuts on vertex position
      if(!((fabs(vertices[0][2]+22.25)<2.25) && (pLead.Mag() < 2.4)))
	continue;

      ep_pmiss_list.push_back(pmiss);

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
      double pmiss_mag = pmiss.Mag();

      // Form the rotated basis
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = oop.Cross(lon).Unit();

      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);

      hist_ep_events->Fill(pmiss_mag,pmiss.Angle(q)*180./M_PI);

      // Now address recoils
      TVector3 prec(pPVec[1][0],pPVec[1][1], pPVec[1][2]);

      // Make Or's cut on recoil vertex position
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > min_prec)))
	{
	  cout << "epp fail\n";
	  continue;
	}

      // Apply Fiducial cuts on the recoil proton
      if (!accept_proton(prec))
	continue;

      TVector3 pcm = prec + pmiss;

      double pcm_lon = pcm.Dot(lon);
      double pcm_inp = pcm.Dot(inp);
      double pcm_oop = pcm.Dot(oop);

      epp_pmiss_list.push_back(pmiss_mag);
      epp_pcm_lon_list.push_back(pcm_lon);
      epp_pcm_inp_list.push_back(pcm_inp);
      epp_pcm_oop_list.push_back(pcm_oop);

      hist_epp_cm_lon->Fill(pmiss_mag,pcm_lon);
      hist_epp_cm_inp->Fill(pmiss_mag,pcm_inp);
      hist_epp_cm_oop->Fill(pmiss_mag,pcm_oop);
      hist_epp_events->Fill(pmiss_mag,pmiss.Angle(q)*180./M_PI);
      hist_epp_recoils->Fill(prec.Mag(),prec.Angle(q)*180./M_PI);
    }
  fepp->Close();

  n_ep_events = ep_pmiss_list.size();
  n_epp_events = epp_pmiss_list.size();
  cerr << "Done reading input data.\n";
  cerr << "Read in " << n_ep_events << " ep events.\n";
  cerr << "Read in " << n_epp_events << " epp events.\n";

  TH1D * hist_ep_pmiss = hist_ep_events->ProjectionX("ep_pmiss");
  int total_samples=0;
  for (int i=0 ; i < n_slices ; i++)
    {
      n_ep_events_per_slice[i] = hist_ep_pmiss->GetBinContent(i+1);
      n_pseudo_per_ep[i] = samples_per_slice / n_ep_events_per_slice[i];
      total_samples += n_pseudo_per_ep[i] * n_ep_events_per_slice[i];
    }

  // Load the acceptance maps
  myMap = new AccMap(argv[3]);

  // Create new memory for parameter guesses
  new_params = new double[5];
  current_params = new double[5];
  random_gauss = new double[3 * total_samples];
  for (int i=0 ; i<3*total_samples ; i++)
    random_gauss[i] = prand->Gaus();

  // Create output rootfile and output tree
  double log_current_post;
  outfile->cd();
  TTree * outtree = new TTree("mcmc","MCMC Samples");
  outtree->Branch("a1",current_params,"a1/D");
  outtree->Branch("a2",current_params+1,"a2/D");
  outtree->Branch("b1",current_params+2,"b1/D");
  outtree->Branch("b2",current_params+3,"b2/D");
  outtree->Branch("sigPerp",current_params+4,"sigPerp/D");
  outtree->Branch("logposterior",&log_current_post,"logposterior/D");

  cerr << "Initializing parameters...\n";
  long double current_post;
  this_sample=-1000;
  do
    {
      initialize_params();
      current_post = posterior(current_params);
      cout << "Attempted initialization: posterior = " << current_post << "\n";
    } while((current_post <= 0) || (isnan(current_post) == true));

  // Loop through MCMC
  int count_true =0;
  for (this_sample=0; this_sample<samples; this_sample++)
    {
      if (this_sample%5 == 0)
	cerr << "Working on iteration " << this_sample << endl;

      pick_new_step();
      long double new_post = posterior(new_params);

      long double post_ratio = new_post/current_post;
      long double proposal_ratio = 1.;

      if (post_ratio*proposal_ratio >= prand->Rndm())
	{
	  // Accept this new step
	  count_true++;
	  current_post = new_post;
	  for (int i = 0; i<5; i++)
	    current_params[i] = new_params[i];
	}

      log_current_post = log(current_post);
      outtree->Fill();
    }

  cerr << "Accept % = " << 100.*count_true/samples << "\n";
  TVectorT <double> acc_pct(1);
  acc_pct[0] = 100.*((double)count_true)/((double)samples);
  acc_pct.Write("acc_pct");

  // Write out histograms
  outfile->cd();
  hist_ep_events->Write();
  hist_epp_cm_lon->Write();
  hist_epp_cm_inp->Write();
  hist_epp_cm_oop->Write();
  hist_epp_events->Write();
  hist_epp_recoils->Write();
  bestLonModel->Write();
  bestInpModel->Write();
  bestOopModel->Write();
  outtree->Write();
  outfile->Close();

  // Clean up
  delete myMap;
  delete prand;
  delete[] new_params;
  delete[] current_params;
  delete[] random_gauss;
  return 0;
}

long double posterior(double * param_set)
{     
  // Safe guard in case a1 or a2 are incorrect
  if((-0.3*param_set[0] + param_set[1]<=0.)||(0.4*param_set[0]+param_set[1]<=0.))
    return 0.;

  // Make pseudo-data set
  TH2D pcm_lon_hist("lon","Longitudinal",n_slices,0.3,1.,n_bins,-pcm_bound,pcm_bound);
  TH2D pcm_inp_hist("inp","Tr. in-plane",n_slices,0.3,1.,n_bins,-pcm_bound,pcm_bound);
  TH2D pcm_oop_hist("oop","Out-of-plane",n_slices,0.3,1.,n_bins,-pcm_bound,pcm_bound);
  TH1D total_acc("acc","Total Accepted",n_slices,0.3,1.);

  int sample_number=0;
  for (int i=0 ; i<n_ep_events ; i++)
    {
      // Establish sigma and mu for the pcm distributions
      double pmiss = ep_pmiss_list[i].Mag();
      double sig_par = param_set[0] * (pmiss - 0.6) + param_set[1];
      double mu_par = param_set[2] * (pmiss - 0.6) + param_set[3];
      double mu_perp = 0.;
      double sig_perp = param_set[4];

      int slice = pcm_lon_hist.GetXaxis()->FindBin(pmiss) - 1;

      for (int j=0 ; j<n_pseudo_per_ep[slice] ; j++, sample_number++)
	{
	  double pcm_lon = mu_par + sig_par*random_gauss[3*sample_number + 0];
	  double pcm_inp = mu_perp+sig_perp*random_gauss[3*sample_number + 1];
	  double pcm_oop = mu_perp+sig_perp*random_gauss[3*sample_number + 2];

	  TVector3 pcm = pcm_lon * ep_lon_list[i] + pcm_inp * ep_inp_list[i] + pcm_oop * ep_oop_list[i];
	  TVector3 prec = pcm - ep_pmiss_list[i];
	  
	  // Fill the histograms with the appropriate acceptance weights
	  double weight = myMap->recoil_accept(prec)/((double)n_pseudo_per_ep[slice]);       
	  pcm_lon_hist.Fill(pmiss,pcm_lon,weight);
	  pcm_inp_hist.Fill(pmiss,pcm_inp,weight);
	  pcm_oop_hist.Fill(pmiss,pcm_oop,weight);
	  total_acc.Fill(pmiss,weight);
	}
    }

  const double scale = ((double) n_epp_events)/((double) n_ep_events);

  pcm_lon_hist.Scale(scale);
  pcm_inp_hist.Scale(scale);
  pcm_oop_hist.Scale(scale);
  total_acc.Scale(scale);

  // pcm_lon_hist.Write();
  //pcm_inp_hist.Write();
  //pcm_oop_hist.Write();
  
  long double result = 1.;
  for (int i=0 ; i<n_epp_events ; i++)
    {
      // Get the measured values of pcm
      double pmiss = epp_pmiss_list[i];
      double pcm_lon = epp_pcm_lon_list[i];
      double pcm_inp = epp_pcm_inp_list[i];
      double pcm_oop = epp_pcm_oop_list[i];
      
      // Get the normalization coefficient
      double norm=total_acc.GetBinContent(total_acc.FindBin(pmiss));

      // Get the model probabilities of pcm
      double pcm_lon_prob = pcm_lon_hist.GetBinContent(pcm_lon_hist.FindBin(pmiss,pcm_lon))/norm;
      double pcm_inp_prob = pcm_inp_hist.GetBinContent(pcm_inp_hist.FindBin(pmiss,pcm_inp))/norm;
      double pcm_oop_prob = pcm_oop_hist.GetBinContent(pcm_oop_hist.FindBin(pmiss,pcm_oop))/norm;

      // Update the result
      if (fabs(pcm_lon) < pcm_bound)
	result *= pcm_lon_prob;
      
      if (fabs(pcm_inp) < pcm_bound)
	result *= pcm_inp_prob;
      
      if (fabs(pcm_oop) < pcm_bound)
	result *= pcm_oop_prob;
    }
  
  // Test if this one is the best so far
  if (result > bestPosterior)
    {
      bestPosterior=result;
      if (bestLonModel)
	{
	  delete bestLonModel;
	  delete bestInpModel;
	  delete bestOopModel;
	}
      char temp_name[100];
      sprintf(temp_name,"lon_%d",this_sample);
      bestLonModel = (TH2D*)pcm_lon_hist.Clone(temp_name);
      sprintf(temp_name,"inp_%d",this_sample);
      bestInpModel = (TH2D*)pcm_inp_hist.Clone(temp_name);
      sprintf(temp_name,"oop_%d",this_sample);
      bestOopModel = (TH2D*)pcm_oop_hist.Clone(temp_name);
    }

  return result;
}

void pick_new_step()
{
  new_params[0] = prand->Gaus(current_params[0], (max_a1-min_a1)/jumpscale);
  new_params[1] = prand->Gaus(current_params[1], (max_a2-min_a2)/jumpscale);
  new_params[2] = prand->Gaus(current_params[2], (max_b1-min_b1)/jumpscale);
  new_params[3] = prand->Gaus(current_params[3], (max_b2-min_b2)/jumpscale);
  new_params[4] = prand->Gaus(current_params[4], (max_sigPerp - min_sigPerp)/jumpscale);
}

void initialize_params()
{
  current_params[0] = prand->Uniform(min_a1,max_a1);
  current_params[1] = prand->Uniform(min_a2,max_a2);
  current_params[2] = prand->Uniform(min_b1,max_b1);
  current_params[3] = prand->Uniform(min_b2,max_b2);
  current_params[4] = prand->Uniform(min_sigPerp,max_sigPerp);

  // First make sure that a1 and a2 are reasonable values
  while((-0.3*current_params[0] + current_params[1]<=0.)||(0.4*current_params[0]+current_params[1]<=0.))
    {
      current_params[0] = prand->Uniform(min_a1,max_a1);
      current_params[1] = prand->Uniform(min_a2,max_a2);
    }
}
 

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

using namespace std;

// Useful math
double sq(double x){return x*x;};
double logGaus(double x, double s, double m);

// Global variables for the parameters
int n_epp_events; 
int n_ep_events;
double * current_params;
double * new_params;
int n_trials = 200.;
double * random_gauss;

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

// Acceptance histograms
TH3D * genHist;
TH3D * accHist;

// Counters for acceptance ratio. To see if our hyper parameters are tuned finely enough
double count_true = 0;
double count_false = 0;
// Random number generator
TRandom3 * prand;

double jumpscale;
long double current_post;

//Starting values for a1, a2, etc
const double a1_bound = .4;
const double a2_bound = .8;
const double b1_bound = .2;
const double b2_bound = .2;
const double sigPerp_bound = .2;
const double pcm_bound = 0.6;

// Helper functions
long double posterior(double * param_set);
void pick_new_step();
bool accept_new_step();
double acceptance_map(TVector3 p);
inline double acceptance(TVector3 p);
double acceptance_fake(TVector3 p);
void initialize_params();


int main(int argc, char ** argv)
{
  if (argc != 6)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "/path/to/data/file /path/to/output/file /path/to/maps.root\n\t[# iter] [jump-scale]\n\n";
      return -1;
    }
  cerr << "Loading data file " << argv[1] << " ...\n";

  // Read in the important info
  int samples = atoi(argv[4]);
  jumpscale = atoi(argv[5]);

  // Initialize random guess
  prand = new TRandom3(0);

  // Read in e'p and e'pp events
  TFile * f = new TFile(argv[1]);
  TTree * t = (TTree*)f->Get("T_all_cuts");
  double qx,qy,qz,p_miss_x,p_miss_y,p_miss_z,p_cm_lon, p_cm_oop, p_cm_inp;
  int nr;
  t->SetBranchAddress("qx",&qx);
  t->SetBranchAddress("qy",&qy);
  t->SetBranchAddress("qz",&qz);
  t->SetBranchAddress("px_p_miss",&p_miss_x);
  t->SetBranchAddress("py_p_miss",&p_miss_y);
  t->SetBranchAddress("pz_p_miss",&p_miss_z);
  t->SetBranchAddress("poofp_p_cm", &p_cm_oop);
  t->SetBranchAddress("ptrans_p_cm",&p_cm_inp);
  t->SetBranchAddress("plong_p_cm", &p_cm_lon);
  t->SetBranchAddress("n_p_r",&nr);
  for (int i=0 ; i<t->GetEntries() ; i++)
    {
      t->GetEvent(i);
      
      TVector3 pmiss(p_miss_x,p_miss_y,p_miss_z);
      ep_pmiss_list.push_back(pmiss);

      // Form the rotated basis
      TVector3 q(qx,qy,qz);
      TVector3 lon = pmiss.Unit();
      TVector3 oop = q.Cross(lon).Unit();
      TVector3 inp = lon.Cross(oop).Unit();
      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);

      // Insist on e'pp events
      if (nr == 1)
        {
	  epp_pmiss_list.push_back(pmiss.Mag());
	  epp_pcm_lon_list.push_back(p_cm_lon);
	  epp_pcm_oop_list.push_back(p_cm_oop);
	  epp_pcm_inp_list.push_back(p_cm_inp);

	}
    }
  n_ep_events = ep_pmiss_list.size();
  n_epp_events = epp_pmiss_list.size();
  f->Close();
  cerr << "Done reading input data.\n";

  // Create new memory for parameter guesses
  new_params = new double[5];
  current_params = new double[5];
  random_gauss = new double[n_ep_events * 3 * n_trials];
  for (int i=0 ; i<n_ep_events*3*n_trials ; i++)
    random_gauss[i] = prand->Gaus();

  // Load the acceptance maps
  TFile * mapFile = new TFile(argv[3]);
  genHist = (TH3D*)mapFile->Get("proton_gen_3D");
  accHist = (TH3D*)mapFile->Get("proton_acc_3D");

  double log_current_post;
  // Create output rootfile and output tree
  cerr << "Creating output file " << argv[2] << " ...\n\n";
  TFile * outfile = new TFile(argv[2],"RECREATE");
  TTree * outtree = new TTree("mcmc","MCMC Samples");
  outtree->Branch("a1",current_params,"a1/D");
  outtree->Branch("a2",current_params+1,"a2/D");
  outtree->Branch("b1",current_params+2,"b1/D");
  outtree->Branch("b2",current_params+3,"b2/D");
  outtree->Branch("sigPerp",current_params+4,"sigPerp/D");
  outtree->Branch("logposterior",&log_current_post,"logposterior/D");

  do
    {
      initialize_params();
      current_post = posterior(current_params);
    } while((current_post <= 0) || (isnan(current_post) == true));

  // Loop through MCMC
  for (int i=0; i<samples; i++)
    {
      if (i%5 == 0)
	cerr << "Working on iteration " << i << endl;

      pick_new_step();

      if (accept_new_step() == true)
	{
	  for (int i = 0; i<5; i++)
	    {
	      current_params[i] = new_params[i];
	    }
	}
      log_current_post = log(current_post);
      outtree->Fill();
    }
  outfile->Write();

  //cout << "Acceptance percentage for model jumpscale of " << jumpscale << " is " << (count_true/samples)*100 << endl;
  TVectorT <double> acceptance_percentage(1);
  acceptance_percentage[0] = (count_true/samples)*100;
  acceptance_percentage.Write("acceptance_percentage");
  delete outfile;
  delete mapFile;
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
  int nBins=20;
  TH2D pcm_lon_hist("lon","Longitudinal",4,0.3,0.7,nBins,-pcm_bound,pcm_bound);
  TH2D pcm_inp_hist("inp","Tr. in-plane",4,0.3,0.7,nBins,-pcm_bound,pcm_bound);
  TH2D pcm_oop_hist("oop","Out-of-plane",4,0.3,0.7,nBins,-pcm_bound,pcm_bound);
  TH1D pcm_lon_hist1d("lon1d","Longitudinal",nBins,-pcm_bound,pcm_bound);
  TH1D pcm_inp_hist1d("inp1d","Tr. in-plane",nBins,-pcm_bound,pcm_bound);
  TH1D pcm_oop_hist1d("oop1d","Out-of-plane",nBins,-pcm_bound,pcm_bound);
  for (int i=0 ; i<n_ep_events ; i++)
    {
      // Establish sigma and mu for the pcm distributions
      double pmiss = ep_pmiss_list[i].Mag();
      double sig_par = param_set[0] * (pmiss - 0.6) + param_set[1];
      double mu_par = param_set[2] * (pmiss - 0.6) + param_set[3];
      double mu_perp = 0.;
      double sig_perp = param_set[4];

      for (int j=0 ; j<n_trials ; j++)
	{
	  double pcm_lon = mu_par + sig_par*random_gauss[3*(n_trials*i + j) + 0];
	  double pcm_inp = mu_perp+sig_perp*random_gauss[3*(n_trials*i + j) + 1];
	  double pcm_oop = mu_perp+sig_perp*random_gauss[3*(n_trials*i + j) + 2];

	  TVector3 pcm = pcm_lon * ep_lon_list[i] + pcm_inp * ep_inp_list[i] + pcm_oop * ep_oop_list[i];
	  TVector3 prec = pcm - ep_pmiss_list[i];
	  
	  // Fill the histograms with the appropriate acceptance weights
	  if (pmiss < 0.7)
	    {
	      pcm_lon_hist.Fill(pmiss,pcm_lon,acceptance(prec));
	      pcm_inp_hist.Fill(pmiss,pcm_inp,acceptance(prec));
	      pcm_oop_hist.Fill(pmiss,pcm_oop,acceptance(prec));
	    }
	  else
	    {
	      pcm_lon_hist1d.Fill(pcm_lon,acceptance(prec));
	      pcm_inp_hist1d.Fill(pcm_inp,acceptance(prec));
	      pcm_oop_hist1d.Fill(pcm_oop,acceptance(prec));
	    }
	}
    }

  const double scale = ((double) n_epp_events)/(((double) n_trials)*((double)n_ep_events));

  pcm_lon_hist.Scale(scale);
  pcm_inp_hist.Scale(scale);
  pcm_oop_hist.Scale(scale);
  pcm_lon_hist1d.Scale(scale);
  pcm_inp_hist1d.Scale(scale);
  pcm_oop_hist1d.Scale(scale);

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
      
      // Get the model probabilities of pcm
      double pcm_lon_model, pcm_inp_model, pcm_oop_model;
      pcm_lon_model = (pmiss < 0.7)? pcm_lon_hist.GetBinContent(pcm_lon_hist.FindBin(pmiss,pcm_lon)) : 
	pcm_lon_hist1d.GetBinContent(pcm_lon_hist.FindBin(pcm_lon));
      pcm_inp_model = (pmiss < 0.7)? pcm_inp_hist.GetBinContent(pcm_inp_hist.FindBin(pmiss,pcm_inp)) : 
	pcm_inp_hist1d.GetBinContent(pcm_inp_hist.FindBin(pcm_inp));
      pcm_oop_model = (pmiss < 0.7)? pcm_oop_hist.GetBinContent(pcm_oop_hist.FindBin(pmiss,pcm_oop)) : 
	pcm_oop_hist1d.GetBinContent(pcm_oop_hist.FindBin(pcm_oop));

      // Update the result
      if (fabs(pcm_lon) < pcm_bound)
	result *= pcm_lon_model;
      
      if (fabs(pcm_inp) < pcm_bound)
	result *= pcm_inp_model;
      
      if (fabs(pcm_oop) < pcm_bound)
	result *= pcm_oop_model;
    }
  return result;
}

void pick_new_step()
{
  new_params[0] = prand->Gaus(current_params[0], a1_bound/jumpscale);
  new_params[1] = prand->Gaus(current_params[1], a2_bound/jumpscale);
  new_params[2] = prand->Gaus(current_params[2], b1_bound/jumpscale);
  new_params[3] = prand->Gaus(current_params[3], b2_bound/jumpscale);
  new_params[4] = prand->Gaus(current_params[4], sigPerp_bound/jumpscale);
}

bool accept_new_step()
{
  current_post = posterior(current_params);
  long double new_post = posterior(new_params);
  long double post_ratio = new_post/current_post;
  //cout << current_post << " -> " << new_post << "    ratio: "<< post_ratio << endl;

  long double proposal_ratio = 1.;
  long double acceptance = post_ratio*proposal_ratio;
  //cout << "Acceptance is: " << acceptance << "\n\n";

  double acceptance_chance = prand->Uniform(1);
  if (acceptance >= acceptance_chance)
    {
      //cout << "True" << endl;
      count_true++;
      current_post = new_post;
      return true;
    }
  else
    {
      //cout << "False" << endl;
      count_false++;
      return false;
    }
}

double resMom(double p)
{
  //The specific parameter values denoted here are obtained by fitting Barak's momentum resolution curve.
  return (.312273 + .287917*p + .550425/(p-.093406))/100;
}

double acceptance_map(TVector3 p)
{
  double mom = p.Mag();
  double cosTheta = p.CosTheta();
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) 
    phi += 360.;

  int momBin = genHist->GetXaxis()->FindBin(mom);
  int cosBin = genHist->GetYaxis()->FindBin(cosTheta);
  int phiBin = genHist->GetZaxis()->FindBin(phi);

  double gen = genHist->GetBinContent(momBin,cosBin,phiBin);
  double acc = accHist->GetBinContent(momBin,cosBin,phiBin);

  if (gen <= 0.)
    return 1.;

  return acc/gen;
}

double acceptance_fake(TVector3 p)
{
  double phi = p.Phi();
  if (phi < -M_PI/6.)
    phi += 2.*M_PI;
  
  int sector = int(phi/60.);
  double san_phi = phi - 60.*sector;
  
  if (fabs(san_phi)<15.)
    return 1.;
  else 
    return 0.;
}

inline double acceptance(TVector3 p)
{
  //return acceptance_map(p);
  return 1;
  //return acceptance_fake(p);
}

void initialize_params()
{
  current_params[0] = prand->Uniform(a1_bound);
  current_params[1] = prand->Uniform(a2_bound);
  current_params[2] = prand->Uniform(b1_bound);
  current_params[3] = prand->Uniform(b2_bound);
  current_params[4] = prand->Uniform(sigPerp_bound);

  // First make sure that a1 and a2 are reasonable values
  while((-0.3*current_params[0] + current_params[1]<=0.)||(0.4*current_params[0]+current_params[1]<=0.))
    {
      current_params[0] = prand->Uniform(a1_bound);
      current_params[1] = prand->Uniform(a2_bound);
    }
}
 

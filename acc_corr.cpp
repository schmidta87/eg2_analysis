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

using namespace std;

// Useful math
double sq(double x){return x*x;};
double logGaus(double x, double s, double m);

// Global variables for the parameters
int n_epp_events; 
int n_ep_events;
int this_sample;
double * current_params;
double * new_params;
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

// Random number generator
TRandom3 * prand;

double jumpscale;

//Starting values for a1, a2, etc
const double a1_bound = .4;
const double a2_bound = .8;
const double b1_bound = .4;
const double b2_bound = .4;
const double sigPerp_bound = .4;
const double pcm_bound = 1.;

// Helper functions
long double posterior(double * param_set);
void pick_new_step();
double acceptance_map(TVector3 p);
inline double acceptance(TVector3 p);
double acceptance_fake(TVector3 p);
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

  // Initialize random guess
  // int seed=0;
  int seed = 777766;
  prand = new TRandom3(seed);
  //prand = new TRandom3(0);
  //seed = prand->Rndm()*1000000;
  //prand->SetSeed(seed);
  cout << "Chosen seed is " << seed << "\n";

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
      TVector3 oop = lon.Cross(TVector3(0.,0.,1.)).Unit();
      TVector3 inp = lon.Cross(oop).Unit();
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
      TVector3 oop = lon.Cross(TVector3(0.,0.,1.)).Unit();
      TVector3 inp = lon.Cross(oop).Unit();

      ep_lon_list.push_back(lon);
      ep_inp_list.push_back(inp);
      ep_oop_list.push_back(oop);

      hist_ep_events->Fill(pmiss_mag,pmiss.Angle(q)*180./M_PI);

      // Now address recoils
      TVector3 prec(pPVec[1][0],pPVec[1][1], pPVec[1][2]);

      // Make Or's cut on recoil vertex position
      if (!((fabs(vertices[1][2]+22.25)<2.25) && (prec.Mag() > 0.35)))
	{
	  cout << "epp fail\n";
	  continue;
	}

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

  // Load the acceptance maps
  TFile * mapFile = new TFile(argv[3]);
  genHist = (TH3D*)mapFile->Get("solid_p_gen");
  accHist = (TH3D*)mapFile->Get("solid_p_acc");
  if ((!genHist)||(!accHist))
    {
      cerr << "The acceptance maps were not loaded properly. Check the histogram names.\n";
      return -2;
    }

  // Create new memory for parameter guesses
  new_params = new double[5];
  current_params = new double[5];
  random_gauss = new double[3 * n_ep_events * n_pseudo];
  for (int i=0 ; i<3*n_ep_events*n_pseudo ; i++)
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
  outtree->Write();
  outfile->Close();

  // Clean up
  mapFile->Close();
  delete prand;
  delete[] new_params;
  delete[] current_params;
  delete[] random_gauss;
  return 0;
}

double logGaus(double x, double m, double s)
{
  return -0.5 * sq((x-m)/s);
}

long double posterior(double * param_set)
{     
  // Safe guard in case a1 or a2 are incorrect
  if((-0.3*param_set[0] + param_set[1]<=0.)||(0.4*param_set[0]+param_set[1]<=0.))
    return 0.;

  vector<double> pcm_lon_list;
  vector<double> pcm_inp_list;
  vector<double> pcm_oop_list;
  vector<double> pmiss_list;
  vector<double> weight_list;

  for (int i=0 ; i<n_ep_events ; i++)
    {
      // Establish sigma and mu for the pcm distributions
      const double pmiss = ep_pmiss_list[i].Mag();
      const double sig_par = param_set[0] * (pmiss - 0.6) + param_set[1];
      const double mu_par = param_set[2] * (pmiss - 0.6) + param_set[3];
      const double mu_perp = 0.;
      const double sig_perp = param_set[4];

      for (int j=0 ; j<n_pseudo ; j++)
	{
	  double pcm_lon = mu_par + sig_par*random_gauss[3*n_pseudo*i + 3*j + 0];
	  double pcm_inp = mu_perp+sig_perp*random_gauss[3*n_pseudo*i + 3*j + 1];
	  double pcm_oop = mu_perp+sig_perp*random_gauss[3*n_pseudo*i + 3*j + 2];

	  TVector3 pcm = pcm_lon * ep_lon_list[i] + pcm_inp * ep_inp_list[i] + pcm_oop * ep_oop_list[i];
	  TVector3 prec = pcm - ep_pmiss_list[i];
	  
	  // Calculate the weight, makin sure the recoils have sufficient momentum
	  double weight;
	  if (prec.Mag()<0.35)
	    weight=0;
	  else
	    weight = acceptance(prec);

	  // Store the pseudo-events in the lists
          pcm_lon_list.push_back(pcm_lon);
          pcm_inp_list.push_back(pcm_inp);
          pcm_oop_list.push_back(pcm_oop);
          pmiss_list.push_back(pmiss);
          weight_list.push_back(weight);	  
	}
    }

  long double sums[n_epp_events];

  #pragma omp parallel for
  for (int i=0 ; i<n_epp_events ; i++)
    {
      sums[i] = 0.;
      for (int j=0; j<n_pseudo*n_ep_events; j++)
        {
	  // Don't waste time on the exponential if out of acceptance
	  if (weight_list[j] <= 0.)
	    continue;

          //Gaussian KDE. No cutoff.
          sums[i] += weight_list[j]*(exp(logGaus(epp_pmiss_list[i],pmiss_list[j],kde_width) + 
					 logGaus(epp_pcm_lon_list[i],pcm_lon_list[j],kde_width) + 
					 logGaus(epp_pcm_inp_list[i],pcm_inp_list[j],kde_width) + 
					 logGaus(epp_pcm_oop_list[i],pcm_oop_list[j],kde_width)));
	}
    }

  // Assemble the result from the product of the sums
  long double result = 1.;
  for (int i=0 ; i<n_epp_events ; i++)
    result *= sums[i];

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

  double gen=0.;
  double acc=0.;
  int exp_size=-1;

  while ((gen <= 0.)&&(exp_size<3))
    {
      exp_size++;
      gen = acc = 0.;
      for (int mB = momBin-exp_size; mB<=momBin+exp_size ; mB++) 
	for (int cB = cosBin-exp_size; cB<=cosBin+exp_size ; cB++) 
	  for (int pB = phiBin-exp_size; pB<=phiBin+exp_size ; pB++) 
	    {
	      gen += genHist->GetBinContent(mB,cB,pB);
	      acc += accHist->GetBinContent(mB,cB,pB);
	    }
    }

  if (gen>0.)
    return acc/gen;
  else
    return 1.;
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
  //return acceptance_fake(p);
  return 1.;
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

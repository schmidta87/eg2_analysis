#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TVectorT.h"

#include "Nuclear_Info.h"
#include "fiducials.h"
#include "helpers.h"

using namespace std;

const double pmiss_cut=0.4;


const double pmiss_lo=0.5;
const double pmiss_md=0.6;
const double pmiss_hi=0.7;

int main(int argc, char ** argv)
{
	if (argc != 4)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   make_hists /path/to/1p/file /path/to/2p/file /path/to/output/file\n\n";
		exit(-1);
	}

	TFile * f1p = new TFile(argv[1]);
	TFile * f2p = new TFile(argv[2]);
	TFile * fo = new TFile(argv[3],"RECREATE");

	// Let's create a vector of all the histogram pointers so we can loop over them, save hassles
	vector<TH1*> h1p_list;
	vector<TH1*> h2p_list;

	// Create histograms
	TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",36,1.,5.);
	h1p_list.push_back(h1p_QSq);
	TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",40,1.2,2.);
	h1p_list.push_back(h1p_xB );
	TH1D * h1p_Pm =  new TH1D("ep_Pm" ,"ep;pMiss [GeV];Counts",30,0.4,1.0);
	h1p_list.push_back(h1p_Pm );
	TH1D * h1p_Pm_coarse =  new TH1D("ep_Pm_coarse" ,"ep;pMiss [GeV];Counts",9,coarse_bin_edges_new);
	h1p_list.push_back(h1p_Pm_coarse);
	TH1D * h1p_Pmq = new TH1D("ep_Pmq","ep;Theta_Pmq [deg];Counts",40,100.,180.);
	h1p_list.push_back(h1p_Pmq);
	TH1D * h1p_cPmq = new TH1D("ep_cPmq","ep;cos(Theta_Pmq);Counts",40,-1.,0.);
	h1p_list.push_back(h1p_cPmq);
	TH1D * h1p_phi1 = new TH1D("ep_phi1","ep;Phi_1 [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phi1);
	TH1D * h1p_phie = new TH1D("ep_phie","ep;Phi_e [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phie);
	TH1D * h1p_theta1 = new TH1D("ep_theta1","ep;Theta_1 [deg];Counts",60,10.,130.);
	h1p_list.push_back(h1p_theta1);
	TH1D * h1p_thetae = new TH1D("ep_thetae","ep;Theta_e [deg];Counts",60,10.,40.);
	h1p_list.push_back(h1p_thetae);
	TH1D * h1p_mome = new TH1D("ep_mome","ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h1p_list.push_back(h1p_mome);
	TH1D * h1p_mom1 = new TH1D("ep_mom1","ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h1p_list.push_back(h1p_mom1);
        TH1D * h1p_alphaq = new TH1D("ep_alphaq","ep;alphaq;Counts",30,-1.5,-0.5);
        h1p_list.push_back(h1p_alphaq);
        TH1D * h1p_alphaLead = new TH1D("ep_alphaLead","ep;alphaLead;Counts",30,0.1,0.6);
        h1p_list.push_back(h1p_alphaLead);
        TH1D * h1p_alphaM = new TH1D("ep_alphaM","ep;alphaM;Counts",30,1.1,1.7);
	h1p_list.push_back(h1p_alphaM);
	TH1D * h1p_Emiss = new TH1D("ep_Emiss","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss);
	TH1D * h1p_Emiss_lo = new TH1D("ep_Emiss_lo","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_lo);
	TH1D * h1p_Emiss_md = new TH1D("ep_Emiss_md","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_md);
	TH1D * h1p_Emiss_hi1 = new TH1D("ep_Emiss_hi1","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi1);
	TH1D * h1p_Emiss_hi2 = new TH1D("ep_Emiss_hi2","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi2);
	TH1D * h1p_Pmq_lo = new TH1D("ep_Pmq_lo","ep;Theta_Pmq [deg];Counts",20,100.,180.);
	h1p_list.push_back(h1p_Pmq_lo);
	TH1D * h1p_Pmq_md = new TH1D("ep_Pmq_md","ep;Theta_Pmq [deg];Counts",20,100.,180.);
	h1p_list.push_back(h1p_Pmq_md);
	TH1D * h1p_Pmq_hi1 = new TH1D("ep_Pmq_hi1","ep;Theta_Pmq [deg];Counts",20,100.,180.);
	h1p_list.push_back(h1p_Pmq_hi1);
	TH1D * h1p_Pmq_hi2 = new TH1D("ep_Pmq_hi2","ep;Theta_Pmq [deg];Counts",20,100.,180.);
	h1p_list.push_back(h1p_Pmq_hi2);
	TH1D * h1p_Pmzq_lo = new TH1D("ep_Pmzq_lo","ep;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h1p_list.push_back(h1p_Pmzq_lo);
	TH1D * h1p_Pmzq_md = new TH1D("ep_Pmzq_md","ep;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h1p_list.push_back(h1p_Pmzq_md);
	TH1D * h1p_Pmzq_hi1 = new TH1D("ep_Pmzq_hi1","ep;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h1p_list.push_back(h1p_Pmzq_hi1);
	TH1D * h1p_Pmzq_hi2 = new TH1D("ep_Pmzq_hi2","ep;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h1p_list.push_back(h1p_Pmzq_hi2);
	TH1D * h1p_PmTq_lo = new TH1D("ep_PmTq_lo","ep;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmTq_lo);
	TH1D * h1p_PmTq_md = new TH1D("ep_PmTq_md","ep;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmTq_md);
	TH1D * h1p_PmTq_hi1 = new TH1D("ep_PmTq_hi1","ep;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmTq_hi1);
	TH1D * h1p_PmTq_hi2 = new TH1D("ep_PmTq_hi2","ep;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmTq_hi2);
	TH1D * h1p_Emiss_fine = new TH1D("ep_Emiss_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_fine);
	TH1D * h1p_Emiss_lo_fine = new TH1D("ep_Emiss_lo_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_lo_fine);
	TH1D * h1p_Emiss_md_fine = new TH1D("ep_Emiss_md_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_md_fine);
	TH1D * h1p_Emiss_hi1_fine = new TH1D("ep_Emiss_hi1_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi1_fine);
	TH1D * h1p_Emiss_hi2_fine = new TH1D("ep_Emiss_hi2_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_hi2_fine);
	TH2D * h1p_pmiss_Emiss = new TH2D("ep_pmiss_Emiss","ep;pmiss [GeV];Emiss [GeV];Counts",24,0.4,1.0,20,-0.2,0.6);
	h1p_list.push_back(h1p_pmiss_Emiss);
	TH2D * h1p_pmiss_E1 = new TH2D("ep_pmiss_E1","ep;pmiss [GeV];E1 [GeV];Counts",24,0.4,1.0,25,0.5,1.0);
	h1p_list.push_back(h1p_pmiss_E1);
	TH2D * h1p_Emiss_by_sector = new TH2D("ep_Emiss_sec","ep;Electron Sector;Emiss [GeV];Counts",6,-0.5,5.5,40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_by_sector);
	TH2D * h1p_pmiss_epsilon = new TH2D("ep_pmiss_epsilon","pmiss_epsilon;pmiss;epsilon;Counts",24,0.4,1.0,25,0.5,1.0);
	h2p_list.push_back(h1p_pmiss_epsilon);
	TH1D * h2p_QSq = new TH1D("epp_QSq","epp;QSq [GeV^2];Counts",30,1.,4.);
	h2p_list.push_back(h2p_QSq);
	TH1D * h2p_xB =  new TH1D("epp_xB" ,"epp;xB;Counts",26,1.2,2.5);
	h2p_list.push_back(h2p_xB );
	TH1D * h2p_Pm =  new TH1D("epp_Pm" ,"epp;pMiss [GeV];Counts",30,0.4,1.0);
	h2p_list.push_back(h2p_Pm );
	TH1D * h2p_Pm_clas =  new TH1D("epp_Pm_clas" ,"epp;pMiss [GeV];Counts",18,0.4,1.0);
	h2p_list.push_back(h2p_Pm );
	TH1D * h2p_Pm_coarse =  new TH1D("epp_Pm_coarse" ,"epp;pMiss [GeV];Counts",9,coarse_bin_edges_new);
	h2p_list.push_back(h2p_Pm_coarse);
	TH1D * h2p_Pmq = new TH1D("epp_Pmq","epp;Theta_Pmq [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmq);
	TH1D * h2p_cPmq = new TH1D("epp_cPmq","epp;cos(Theta_Pmq);Counts",20,-1.,0.);
	h2p_list.push_back(h2p_cPmq);
	TH1D * h2p_Pmr = new TH1D("epp_Pmr","epp;Theta_Pmr [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmr);
	TH1D * h2p_cPmr = new TH1D("epp_cPmr","epp;cos(Theta_Pmr);Counts",20,-1.,0.);
	h2p_list.push_back(h2p_cPmr);
	TH1D * h2p_phi1 = new TH1D("epp_phi1","epp;Phi_1 [deg];Counts",60,-30.,330.);
	h2p_list.push_back(h2p_phi1);
	TH1D * h2p_phi2 = new TH1D("epp_phi2","epp;Phi_2 [deg];Counts",60,-30.,330.);
	h2p_list.push_back(h2p_phi2);
	TH1D * h2p_phie = new TH1D("epp_phie","epp;Phi_e [deg];Counts",60,-30.,330.);
	h2p_list.push_back(h2p_phie);
	TH1D * h2p_thetae = new TH1D("epp_thetae","epp;Theta_e [deg];Counts",30,10.,40.);
	h2p_list.push_back(h2p_thetae);
	TH1D * h2p_mome = new TH1D("epp_mome","epp;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h2p_list.push_back(h2p_mome);
	TH1D * h2p_theta1 = new TH1D("epp_theta1","epp;Theta_1 [deg];Counts",30,10.,130.);
	h2p_list.push_back(h2p_theta1);
	TH1D * h2p_theta2 = new TH1D("epp_theta2","epp;Theta_2 [deg];Counts",30,10.,130.);
	h2p_list.push_back(h2p_theta2);
	TH1D * h2p_mom2 = new TH1D("epp_mom2","epp;Recoil Mom [GeV/c];Counts",17,0.35,1.2);
	h2p_list.push_back(h2p_mom2);
	TH1D * h2p_mom1 = new TH1D("epp_mom1","epp;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h2p_list.push_back(h2p_mom1);
	TH1D * h2p_alphaq = new TH1D("epp_alphaq","ep;alphaq;Counts",30,-1.5,-0.5);
        h2p_list.push_back(h2p_alphaq);
        TH1D * h2p_alphaLead = new TH1D("epp_alphaLead","ep;alphaLead;Counts",30,0.1,0.6);
        h2p_list.push_back(h2p_alphaLead);
        TH1D * h2p_alphaM = new TH1D("epp_alphaM","ep;alphaM;Counts",30,1.1,1.7);
        h2p_list.push_back(h2p_alphaM);
        TH1D * h2p_alphaRec = new TH1D("epp_alphaRec","ep;alphaRec;Counts",30,0.4,1.3);
        h2p_list.push_back(h2p_alphaRec);
        TH1D * h2p_alphaD = new TH1D("epp_alphaD","ep;alphaD;Counts",30,1.6,2.6);
        h2p_list.push_back(h2p_alphaD);
	TH1D * h2p_Emiss = new TH1D("epp_Emiss","epp;Emiss [GeV];Counts",40,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss);
	TH1D * h2p_Emiss_lo = new TH1D("epp_Emiss_lo","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_lo);
	TH1D * h2p_Emiss_md = new TH1D("epp_Emiss_md","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_md);
	TH1D * h2p_Emiss_hi1 = new TH1D("epp_Emiss_hi1","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi1);
	TH1D * h2p_Emiss_hi2 = new TH1D("epp_Emiss_hi2","epp;Emiss [GeV];Counts",20,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi2);
	TH1D * h2p_Emiss_fine = new TH1D("epp_Emiss_fine","epp;Emiss [GeV];Counts",160,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_fine);
	TH1D * h2p_Emiss_lo_fine = new TH1D("epp_Emiss_lo_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_lo_fine);
	TH1D * h2p_Emiss_md_fine = new TH1D("epp_Emiss_md_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_md_fine);
	TH1D * h2p_Emiss_hi1_fine = new TH1D("epp_Emiss_hi1_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi1_fine);
	TH1D * h2p_Emiss_hi2_fine = new TH1D("epp_Emiss_hi2_fine","epp;Emiss [GeV];Counts",80,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_hi2_fine);
	TH1D * h2p_Pmq_lo = new TH1D("epp_Pmq_lo","epp;Theta_Pmq [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmq_lo);
	TH1D * h2p_Pmq_md = new TH1D("epp_Pmq_md","epp;Theta_Pmq [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmq_md);
	TH1D * h2p_Pmq_hi1 = new TH1D("epp_Pmq_hi1","epp;Theta_Pmq [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmq_hi1);
	TH1D * h2p_Pmq_hi2 = new TH1D("epp_Pmq_hi2","epp;Theta_Pmq [deg];Counts",20,100.,180.);
	h2p_list.push_back(h2p_Pmq_hi2);
	TH1D * h2p_Pmzq_lo = new TH1D("epp_Pmzq_lo","epp;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h2p_list.push_back(h2p_Pmzq_lo);
	TH1D * h2p_Pmzq_md = new TH1D("epp_Pmzq_md","epp;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h2p_list.push_back(h2p_Pmzq_md);
	TH1D * h2p_Pmzq_hi1 = new TH1D("epp_Pmzq_hi1","epp;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h2p_list.push_back(h2p_Pmzq_hi1);
	TH1D * h2p_Pmzq_hi2 = new TH1D("epp_Pmzq_hi2","epp;pMiss_||q [GeV];Counts",30,-1.0,0.0);
	h2p_list.push_back(h2p_Pmzq_hi2);
	TH1D * h2p_PmTq_lo = new TH1D("epp_PmTq_lo","epp;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h2p_list.push_back(h2p_PmTq_lo);
	TH1D * h2p_PmTq_md = new TH1D("epp_PmTq_md","epp;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h2p_list.push_back(h2p_PmTq_md);
	TH1D * h2p_PmTq_hi1 = new TH1D("epp_PmTq_hi1","epp;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h2p_list.push_back(h2p_PmTq_hi1);
	TH1D * h2p_PmTq_hi2 = new TH1D("epp_PmTq_hi2","epp;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h2p_list.push_back(h2p_PmTq_hi2);
	TH2D * h2p_pmiss_E1 = new TH2D("epp_pmiss_E1","epp;pmiss [GeV];E1 [GeV];Counts",24,0.4,1.0,25,0.5,1.0);
	h2p_list.push_back(h2p_pmiss_E1);
	TH2D * h2p_pmiss_appEstar = new TH2D("epp_pmiss_appEstar","epp;pmiss [GeV];Apparent Estar [GeV];Counts",24,0.4,1.0,20,-0.2,0.8);
	h2p_list.push_back(h2p_pmiss_appEstar);
	TH2D * h2p_pmiss_epsilon = new TH2D("epp_pmiss_epsilon","pmiss_epsilon;pmiss;epsilon;Counts",24,0.4,1.0,25,0.5,1.0);
	h2p_list.push_back(h2p_pmiss_epsilon);
	TH2D * h2p_pRec_epsilon = new TH2D("pRec_epsilon","pRec_epsilon;pRec;epsilon;Counts",20,0.3,1.2,20,0.2,1.2);
	h2p_list.push_back(h2p_pRec_epsilon);
	TH2D * h2p_pRec_eMiss = new TH2D("pRec_eMiss","pRec_eMiss;pRec;eMiss;Counts",10,0.35,0.9,10,0,0.5);
	h2p_list.push_back(h2p_pRec_eMiss);

	TH2D * pp_to_p_2d = new TH2D("pp_to_p_2d","2d ratio;pmiss [GeV];E1 [GeV];pp/p",28,0.35,1.0,20,0.5,0.9);

	TH1D * h1p_Emiss_byBin[12];
	TH1D * h2p_Emiss_byBin[12];
	for (int i=0 ; i<12 ; i++)
	  {
	    char temp[100];

	    sprintf(temp,"ep_Emiss_%d",i);
	    h1p_Emiss_byBin[i] = new TH1D(temp,"ep;Emiss [GeV];Counts",160,-0.2,0.6);
	    h1p_list.push_back(h1p_Emiss_byBin[i]);

	    sprintf(temp,"epp_Emiss_%d",i);
	    h2p_Emiss_byBin[i] = new TH1D(temp,"epp;Emiss [GeV];Counts",120,-0.2,0.6);
	    h2p_list.push_back(h2p_Emiss_byBin[i]);
	  }



	TH1D * h1p_thetae_bySec[6];
	TH1D * h1p_theta1_bySec[6];
	TH1D * h2p_theta1_bySec[6];
	TH1D * h2p_theta2_bySec[6];
	for (int i=0 ; i<6 ; i++)
	  {
	    char temp[100];
	    
	    sprintf(temp,"ep_theta1_%d",i);
	    h1p_theta1_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,130.);
	    h1p_list.push_back(h1p_theta1_bySec[i]);
	    
	    sprintf(temp,"ep_thetae_%d",i);
	    h1p_thetae_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,40.);
	    h1p_list.push_back(h1p_thetae_bySec[i]);
	    
	    sprintf(temp,"epp_theta1_%d",i);
	    h2p_theta1_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
	    h2p_list.push_back(h2p_theta1_bySec[i]);
	    
	    sprintf(temp,"epp_theta2_%d",i);
	    h2p_theta2_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
	    h2p_list.push_back(h2p_theta2_bySec[i]);
	  }

	// Now that all histograms have been defined, set them to Sumw2
	for (int i=0 ; i<h1p_list.size() ; i++)
	  h1p_list[i]->Sumw2();
	for (int i=0 ; i<h2p_list.size() ; i++)
	  h2p_list[i]->Sumw2();

	// pp2p graphs
	TGraphAsymmErrors * pp_to_p = new TGraphAsymmErrors();
	pp_to_p->SetName("pp_to_p");
	pp_to_p->SetTitle("pp_to_p;p_miss [GeV];pp_to_p ratio");
	TGraphAsymmErrors * pp_to_p_coarse = new TGraphAsymmErrors();
	pp_to_p_coarse->SetName("pp_to_p_coarse");
	pp_to_p_coarse->SetTitle("pp_to_p;p_miss [GeV];pp_to_p ratio");

	// Loop over 1p tree
	cerr << " Looping over 1p tree...\n";
	TTree * t1p = (TTree*)f1p->Get("T");
	Float_t Xb, Q2, Pmiss_size[2], Pp[2][3], Rp[2][3], Pp_size[2], Pmiss_q_angle[2], Pe[3], q[3];
	Double_t weight = 1.;
	t1p->SetBranchAddress("Pmiss_q_angle",Pmiss_q_angle);
	t1p->SetBranchAddress("Xb",&Xb);
	t1p->SetBranchAddress("Q2",&Q2);
	t1p->SetBranchAddress("Pmiss_size",Pmiss_size);
	t1p->SetBranchAddress("Pp_size",Pp_size);
	t1p->SetBranchAddress("Rp",Rp);
	t1p->SetBranchAddress("Pp",Pp);
	t1p->SetBranchAddress("Pe",Pe);
	t1p->SetBranchAddress("q",q);

	// See if there is a weight branch
	TBranch * weight_branch = t1p->GetBranch("weight");
	if (weight_branch)
	{
		t1p->SetBranchAddress("weight",&weight);
	}

	for (int event =0 ; event < t1p->GetEntries() ; event++)
	{
		t1p->GetEvent(event);

		// Do necessary cuts
		if (fabs(Rp[0][2]+22.25)>2.25)
		  continue;
		if (Pp_size[0]>2.4)
		  continue;
		if (Pmiss_size[0]<pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;
		
                // A few more vectors                                                            
                TVector3 vq(q[0],q[1],q[2]);
                TVector3 vqUnit = vq.Unit();
		TVector3 vm = vp - vq;

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		// Kinematic variables we need
		double omega = Q2/(2.*mN*Xb);
		double Ep = sqrt(Pp_size[0]*Pp_size[0] + mN*mN);
		double Emiss = -m_12C + mN + sqrt( sq(omega + m_12C - Ep) - (Pmiss_size[0]*Pmiss_size[0]));
		double epsilon = Ep - omega;

		//Let's calculate light cone variables                                           
                double alphaq= (omega - vq.Mag()) / mN;
                double alphaLead = (Ep - vp.Dot(vqUnit)) / mN;
                double alphaM = alphaLead - alphaq;

                h1p_alphaq->Fill(alphaq,weight);
                h1p_alphaLead->Fill(alphaLead,weight);
                h1p_alphaM->Fill(alphaM,weight);

		// Let's make a sanitized phi and sector
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = clas_sector(phie_deg);
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome->Fill(ve.Mag(),weight);

		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vp.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		
		// Let's figure out missing energy! 
		h1p_Emiss->Fill(Emiss,weight);
		h1p_Emiss_fine->Fill(Emiss,weight);
		h1p_pmiss_Emiss->Fill(Pmiss_size[0],Emiss,weight);
		h1p_Emiss_by_sector->Fill(sec_e,Emiss);

		if (Pmiss_size[0] < pmiss_lo)
		  {
		    h1p_Emiss_lo->Fill(Emiss,weight);
		    h1p_Emiss_lo_fine->Fill(Emiss,weight);
		    h1p_Pmq_lo->Fill(Pmiss_q_angle[0],weight);
		    h1p_Pmzq_lo->Fill(vm.Dot(vqUnit),weight);
		    h1p_PmTq_lo->Fill(vm.Perp(vqUnit),weight);
		  }
		else if (Pmiss_size[0] < pmiss_md)
		  {
		    h1p_Emiss_md->Fill(Emiss,weight);
		    h1p_Emiss_md_fine->Fill(Emiss,weight);
		    h1p_Pmq_md->Fill(Pmiss_q_angle[0],weight);
		    h1p_Pmzq_md->Fill(vm.Dot(vqUnit),weight);
		    h1p_PmTq_md->Fill(vm.Perp(vqUnit),weight);
		  }
		else if (Pmiss_size[0] < 1.)
		  {
		    if (Pmiss_size[0] < pmiss_hi)
		      {
			h1p_Emiss_hi1->Fill(Emiss,weight);
			h1p_Emiss_hi1_fine->Fill(Emiss,weight);
			h1p_Pmq_hi1->Fill(Pmiss_q_angle[0],weight);
			h1p_Pmzq_hi1->Fill(vm.Dot(vqUnit),weight);
			h1p_PmTq_hi1->Fill(vm.Perp(vqUnit),weight);
		      }
		    else
		      {
			h1p_Emiss_hi2->Fill(Emiss,weight);
			h1p_Emiss_hi2_fine->Fill(Emiss,weight);
			h1p_Pmq_hi2->Fill(Pmiss_q_angle[0],weight);
			h1p_Pmzq_hi2->Fill(vm.Dot(vqUnit),weight);
			h1p_PmTq_hi2->Fill(vm.Perp(vqUnit),weight);
		      }
		  }

		h1p_pmiss_E1->Fill(Pmiss_size[0],epsilon,weight);

		int pmiss_bin = (Pmiss_size[0]-0.4)/0.05;
		h1p_Emiss_byBin[pmiss_bin]->Fill(Emiss,weight);

		h1p_pmiss_epsilon->Fill(Pmiss_size[0],epsilon,weight);
	}

	// Loop over 2p tree
	cerr << " Looping over 2p tree...\n";
	TTree * t2p = (TTree*)f2p->Get("T");
	t2p->SetBranchAddress("Pmiss_q_angle",Pmiss_q_angle);
	t2p->SetBranchAddress("Xb",&Xb);
	t2p->SetBranchAddress("Q2",&Q2);
	t2p->SetBranchAddress("Pmiss_size",Pmiss_size);
	t2p->SetBranchAddress("Pp_size",Pp_size);
	t2p->SetBranchAddress("Rp",Rp);
	t2p->SetBranchAddress("Pp",Pp);
	t2p->SetBranchAddress("Pe",Pe);
	t2p->SetBranchAddress("q",q);
	// See if there is a weight branch
	weight=1.;
	weight_branch = t2p->GetBranch("weight");
	if (weight_branch)
	{
		t2p->SetBranchAddress("weight",&weight);
	}
	for (int event =0 ; event < t2p->GetEntries() ; event++)
	{
		t2p->GetEvent(event);

		// Do necessary cuts
		if (fabs(Rp[0][2]+22.25)>2.25)
			continue;
		if (Pp_size[0]>2.4)
			continue;
		if (Pmiss_size[0]<pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vlead))
		  continue;

                // A few more vectors                                                             
                TVector3 vq(q[0],q[1],q[2]);
                TVector3 vqUnit = vq.Unit();
                TVector3 vmiss = vlead - vq;

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		// Kinematic variables we need
		double omega = Q2/(2.*mN*Xb);
		double Elead = sqrt(Pp_size[0]*Pp_size[0] + mN*mN);
		double Emiss = -m_12C + mN + sqrt( sq(omega + m_12C - Elead) - (Pmiss_size[0]*Pmiss_size[0]));
		double epsilon = Elead - omega;

                //Let's calculate light cone variables                                                      
                double alphaq= (omega - vq.Mag()) / mN;
                double alphaLead = (Elead - vlead.Dot(vqUnit)) / mN;
                double alphaM = alphaLead - alphaq;

                h2p_alphaq->Fill(alphaq,weight);
                h2p_alphaLead->Fill(alphaLead,weight);
                h2p_alphaM->Fill(alphaM,weight);

		// Let's make a sanitized phi and sector
		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = clas_sector(phie_deg);
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);
		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome->Fill(ve.Mag(),weight);

		double phi1_deg = vlead.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vlead.Theta() * 180./M_PI;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		
		h1p_Emiss->Fill(Emiss,weight);
		h1p_Emiss_fine->Fill(Emiss,weight);
		h1p_pmiss_Emiss->Fill(Pmiss_size[0],Emiss,weight);
		h1p_pmiss_E1->Fill(Pmiss_size[0],epsilon,weight);
		h1p_Emiss_by_sector->Fill(sec_e,Emiss);
		if (Pmiss_size[0] < pmiss_lo)
		  {
		    h1p_Emiss_lo->Fill(Emiss,weight);
		    h1p_Emiss_lo_fine->Fill(Emiss,weight);
		    h1p_Pmq_lo->Fill(Pmiss_q_angle[0],weight);
		    h1p_Pmzq_lo->Fill(vmiss.Dot(vqUnit),weight);
		    h1p_PmTq_lo->Fill(vmiss.Perp(vqUnit),weight);
		  }
		else if (Pmiss_size[0] < pmiss_md)
		  {
		    h1p_Emiss_md->Fill(Emiss,weight);
		    h1p_Emiss_md_fine->Fill(Emiss,weight);
		    h1p_Pmq_md->Fill(Pmiss_q_angle[0],weight);
		    h1p_Pmzq_md->Fill(vmiss.Dot(vqUnit),weight);
		    h1p_PmTq_md->Fill(vmiss.Perp(vqUnit),weight);
		  }
		else if (Pmiss_size[0] < 1.)
		  {
		    if (Pmiss_size[0] < pmiss_hi)
		      {
			h1p_Emiss_hi1->Fill(Emiss,weight);
			h1p_Emiss_hi1_fine->Fill(Emiss,weight);
			h1p_Pmq_hi1->Fill(Pmiss_q_angle[0],weight);
			h1p_Pmzq_hi1->Fill(vmiss.Dot(vqUnit),weight);
			h1p_PmTq_hi1->Fill(vmiss.Perp(vqUnit),weight);
		      }
		    else
		      {
			h1p_Emiss_hi2->Fill(Emiss,weight);
			h1p_Emiss_hi2_fine->Fill(Emiss,weight);
			h1p_Pmq_hi2->Fill(Pmiss_q_angle[0],weight);
			h1p_Pmzq_hi2->Fill(vmiss.Dot(vqUnit),weight);
			h1p_PmTq_hi2->Fill(vmiss.Perp(vqUnit),weight);
		      }
		  }

		int pmiss_bin = (Pmiss_size[0]-0.4)/0.05;
		h1p_Emiss_byBin[pmiss_bin]->Fill(Emiss,weight);

		h1p_pmiss_epsilon->Fill(Pmiss_size[0],epsilon,weight);

		// Make a check on the recoils
		if (fabs(Rp[1][2]+22.25)>2.25)
		  continue;
		if (Pp_size[1] < 0.35)
		  continue;
		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);      
		if (!accept_proton(vrec))
		  continue;

                double Erec = sqrt(vrec.Mag2() + mN*mN);
                 TVector3 vcm = vmiss + vrec;

                // Find final light cone variables                                                          
                double alphaRec = (Erec - vrec.Dot(vqUnit))/mN;
                double alphaD = alphaM + alphaRec;

		h2p_QSq->Fill(Q2,weight);
		h2p_xB ->Fill(Xb,weight);
		h2p_Pm ->Fill(Pmiss_size[0],weight);
		h2p_Pm_clas ->Fill(Pmiss_size[0],weight);
		h2p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h2p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h2p_Pmr->Fill(vmiss.Angle(vrec)*180./M_PI,weight);
		h2p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);
		h2p_cPmr->Fill(cos(vmiss.Angle(vrec)),weight);

                h2p_alphaRec->Fill(alphaRec,weight);
                h2p_alphaD->Fill(alphaD,weight);

		h2p_phie->Fill(phie_deg,weight);
		h2p_thetae->Fill(thetae_deg,weight);
		h2p_mome->Fill(ve.Mag(),weight);

		h2p_phi1->Fill(phi1_deg,weight);
		h2p_theta1->Fill(theta1_deg,weight);
		h2p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h2p_mom1->Fill(Pp_size[0],weight);

		h2p_Emiss->Fill(Emiss,weight);
		h2p_Emiss_fine->Fill(Emiss,weight);
		h2p_pmiss_E1->Fill(Pmiss_size[0],epsilon,weight);

		h2p_pmiss_epsilon->Fill(Pmiss_size[0],epsilon,weight);
		h2p_pRec_epsilon->Fill(Pp_size[1],epsilon,weight);
	        h2p_pRec_eMiss->Fill(Pp_size[1],Emiss,weight);

		if (Pmiss_size[0] < pmiss_lo)
		  {
		    h2p_Emiss_lo->Fill(Emiss,weight);
		    h2p_Emiss_lo_fine->Fill(Emiss,weight);
		    h2p_Pmq_lo->Fill(Pmiss_q_angle[0],weight);
		    h2p_Pmzq_lo->Fill(vmiss.Dot(vqUnit),weight);
		    h2p_PmTq_lo->Fill(vmiss.Perp(vqUnit),weight);
		  }
		else if (Pmiss_size[0] < pmiss_md)
		  {
		    h2p_Emiss_md->Fill(Emiss,weight);
		    h2p_Emiss_md_fine->Fill(Emiss,weight);
		    h2p_Pmq_md->Fill(Pmiss_q_angle[0],weight);
		    h2p_Pmzq_md->Fill(vmiss.Dot(vqUnit),weight);
		    h2p_PmTq_md->Fill(vmiss.Perp(vqUnit),weight);
		  }
		else if (Pmiss_size[0] < 1.)
		  {
		    if (Pmiss_size[0] < pmiss_hi)
		      {
			h2p_Emiss_hi1->Fill(Emiss,weight);
			h2p_Emiss_hi1_fine->Fill(Emiss,weight);
			h2p_Pmq_hi1->Fill(Pmiss_q_angle[0],weight);
			h2p_Pmzq_hi1->Fill(vmiss.Dot(vqUnit),weight);
			h2p_PmTq_hi1->Fill(vmiss.Perp(vqUnit),weight);
		      }
		    else
		      {
			h2p_Emiss_hi2->Fill(Emiss,weight);
			h2p_Emiss_hi2_fine->Fill(Emiss,weight);
			h2p_Pmq_hi2->Fill(Pmiss_q_angle[0],weight);
			h2p_Pmzq_hi2->Fill(vmiss.Dot(vqUnit),weight);
			h2p_PmTq_hi2->Fill(vmiss.Perp(vqUnit),weight);
		      }
		  }

		// Let's make a sanitized phi and sector
		double phi2_deg = vrec.Phi() * 180./M_PI;
		if (phi2_deg < -30.)
			phi2_deg += 360.;
		int sector2 = clas_sector(phi2_deg);
		double theta2_deg = vrec.Theta() * 180./M_PI;

		h2p_phi2->Fill(phi2_deg,weight);
		h2p_theta2->Fill(theta2_deg,weight);
		h2p_theta2_bySec[sector2]->Fill(theta2_deg,weight);
		h2p_mom2->Fill(Pp_size[1],weight);

		// Histogram for the "apparent E*"
		double apparent_Estar = sqrt(sq(sqrt(sq(m_10B) + vcm.Mag2()) + Erec) -vlead.Mag2()) - m_11B;
		h2p_pmiss_appEstar->Fill(Pmiss_size[0],apparent_Estar,weight);

		h2p_Emiss_byBin[pmiss_bin]->Fill(Emiss,weight);
	}
	f1p->Close();
	f2p->Close();

	cerr << "The ep and epp integrals are: " << h1p_Pm->Integral() << " "  << h2p_Pm->Integral() << "\n";
	cerr << "Broken down by pmiss range...\n\n";
	for (int j=4 ; j<=24 ; j+=4)
	  {
	    double min=0.4 + 0.1*(j-4)/4.;
	    double max=0.4 + 0.1*j/4.;
	    cerr << min << " < pmiss < " << max << " : " << h1p_Pm->Integral(j-3,j) << " " << h2p_Pm->Integral(j-3,j) << "\n";
	  }

	// pp-to-p
	pp_to_p->BayesDivide(h2p_Pm,h1p_Pm);
	pp_to_p_coarse->BayesDivide(h2p_Pm_coarse,h1p_Pm_coarse);
  	for (int binX=1 ; binX<=pp_to_p_2d->GetNbinsX() ; binX++)
	  for (int binY=1 ; binY<=pp_to_p_2d->GetNbinsY() ; binY++)
	    {
	      
	      double N_pp = h2p_pmiss_E1->GetBinContent(binX,binY);
	      double N_p =  h1p_pmiss_E1->GetBinContent(binX,binY);
	      double ratio = 0.;
	      if (N_p >0)
		ratio = N_pp / N_p;
	      pp_to_p_2d->SetBinContent(binX,binY,ratio);
	    }

	// Write out
	fo -> cd();
	pp_to_p->Write();
	pp_to_p_coarse->Write();
	pp_to_p_2d->Write();
	
	const double data_ep = 5854.;
	const double data_ep_cor = 6077.;
	const double data_epp = 411.;
	const double pnorm = data_ep/h1p_Pm->Integral();
	const double ppnorm = pnorm;

	// Including a factor if we watn to rescale data to match epp luminosity
	const double renorm = data_epp/h2p_Pm->Integral()/(data_ep/h1p_Pm->Integral());
	TVectorT<double> renorm_factor(1);
	renorm_factor[0] = renorm;
	renorm_factor.Write("factor");

	// scale all the histograms, and write them out
	for (int i=0 ; i<h1p_list.size() ; i++)
	  {
	    h1p_list[i]->Scale(pnorm);
	    h1p_list[i]->Write();
	  }
	for (int i=0 ; i<h2p_list.size() ; i++)
	  {
	    h2p_list[i]->Scale(ppnorm);
	    h2p_list[i]->Write();
	  }

	fo->Close();

	return 0;
}

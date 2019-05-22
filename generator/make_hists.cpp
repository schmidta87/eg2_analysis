#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TVectorT.h"

#include "Nuclear_Info.h"
#include "Cross_Sections.h"
#include "fiducials.h"
#include "helpers.h"
#include "AccMap.h"

using namespace std;

const double pmiss_cut=0.4;

const double pmiss_lo=0.5;
const double pmiss_md=0.6;
const double pmiss_hi=0.7;
const double Ebeam=eg2beam;
const double sCutOff=0;
const double sigmaCM=0.143;
const double Estar=0.03;

const double acc_thresh=0.0;

int main(int argc, char ** argv)
{
	if (argc < 5)
	{
		cerr << "Wrong number of arguments. Instead try:\n"
			<< "   make_hists /path/to/1p/file /path/to/2p/file /path/to/map/file /path/to/output/file\n\n";
		exit(-1);
	}
	
	TFile * f1p = new TFile(argv[1]);
	TFile * f2p = new TFile(argv[2]);
	TFile * fo = new TFile(argv[4],"RECREATE");
	Cross_Sections myCS(cc1,kelly);
	bool doSWeight = false;
	bool doCut = true;
	bool doOtherCut = true;
	bool doGaps = true;
	int c;
	while((c=getopt (argc-4, &argv[4], "SCOg")) != -1)
	  switch(c)
	    {
	    case 'S':
	      doSWeight = true;
	      break;
	    case 'C':
	      doCut = false;
	      break;
	    case 'O':
	      doOtherCut = false;
	      break;
	    case 'g':
	      doGaps = false;
	      break;
	    case '?':
	      return -1;
	    default:
	      abort();
	    }

	// We'll need to get acceptance maps in order to do a fiducial cut on minimum acceptance
	AccMap proton_map(argv[3], "p");

	// Create some directories so that the output file isn't so crowded in a browser
	TDirectory * dir_by_sec = fo->mkdir("by_sec");
	TDirectory * dir_sub_bins = fo->mkdir("sub_bins");
	fo->cd();

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
	TH1D * h1p_dE1 = new TH1D("ep_dE1","ep;dE1;Counts",40,-0.2,1);
	h1p_list.push_back(h1p_dE1);
	TH1D * h1p_rE1 = new TH1D("ep_rE1","ep;rE1;Counts",40,0,1);
	h1p_list.push_back(h1p_rE1);		
	TH2D * h1p_pmiss_Emiss = new TH2D("ep_pmiss_Emiss","ep;pmiss [GeV];Emiss [GeV];Counts",24,0.4,1.0,20,-0.2,0.6);
	h1p_list.push_back(h1p_pmiss_Emiss);
	TH2D * h1p_pmiss_E1 = new TH2D("ep_pmiss_E1","ep;pmiss [GeV];E1 [GeV];Counts",24,0.4,1.0,25,0.5,1.0);
	h1p_list.push_back(h1p_pmiss_E1);
	TH2D * h1p_Emiss_by_sector = new TH2D("ep_Emiss_sec","ep;Electron Sector;Emiss [GeV];Counts",6,-0.5,5.5,40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_by_sector);
	TH2D * h1p_pmiss_epsilon = new TH2D("ep_pmiss_epsilon","pmiss_epsilon;pmiss;epsilon;Counts",24,0.4,1.0,25,0.5,1.0);
	h1p_list.push_back(h1p_pmiss_epsilon);
	TH1D * h2p_QSq = new TH1D("epp_QSq","epp;QSq [GeV^2];Counts",30,1.,4.);
	h2p_list.push_back(h2p_QSq);
	TH1D * h2p_xB =  new TH1D("epp_xB" ,"epp;xB;Counts",26,1.2,2.5);
	h2p_list.push_back(h2p_xB );
	TH1D * h2p_Pm =  new TH1D("epp_Pm" ,"epp;pMiss [GeV];Counts",30,0.4,1.0);
	h2p_list.push_back(h2p_Pm );
	TH1D * h2p_Pm_clas =  new TH1D("epp_Pm_clas" ,"epp;pMiss [GeV];Counts",18,0.4,1.0);
	h2p_Pm_clas->Sumw2();
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
	TH1D * h2p_momdiff = new TH1D("epp_momdiff","epp;Delta_mom [GeV/c];Counts",40,0.0,2.0);
	h2p_list.push_back(h2p_momdiff);
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
	TH1D * h2p_dE1 = new TH1D("epp_dE1","epp;dE1;Counts",40,-0.2,1);
	h2p_list.push_back(h2p_dE1);
	TH1D * h2p_rE1 = new TH1D("epp_rE1","epp;rE1;Counts",40,0,1);
	h2p_list.push_back(h2p_rE1);		
	TH2D * h2p_pmiss_E1 = new TH2D("epp_pmiss_E1","epp;pmiss [GeV];E1 [GeV];Counts",24,0.4,1.0,25,0.5,1.0);
	h2p_list.push_back(h2p_pmiss_E1);
	TH2D * h2p_pmiss_appEstar = new TH2D("epp_pmiss_appEstar","epp;pmiss [GeV];Apparent Estar [GeV];Counts",24,0.4,1.0,20,-0.2,0.8);
	h2p_list.push_back(h2p_pmiss_appEstar);
	TH2D * h2p_pmiss_epsilon = new TH2D("epp_pmiss_epsilon","pmiss_epsilon;pmiss;epsilon;Counts",24,0.4,1.0,25,0.5,1.0);
	h2p_list.push_back(h2p_pmiss_epsilon);

	TH2D * h1p_pmiss_QSq = new TH2D("ep_pmiss_QSq","pmiss_epsilon;pmiss;Q^2;Counts",24,0.4,1.0,30,1.,4.);
	h1p_list.push_back(h1p_pmiss_QSq);
	TH2D * h2p_pmiss_QSq = new TH2D("epp_pmiss_QSq","pmiss_epsilon;pmiss;Q^2;Counts",24,0.4,1.0,30,1.,4.);
	h2p_list.push_back(h2p_pmiss_QSq);
	TH2D * h1p_pmiss_mom1 = new TH2D("ep_pmiss_mom1","pmiss_epsilon;pmiss;pLead;Counts",24,0.4,1.0,40,0.4,2.4);
	h1p_list.push_back(h1p_pmiss_mom1);
	TH2D * h2p_pmiss_mom1 = new TH2D("epp_pmiss_mom1","pmiss_epsilon;pmiss;pLead;Counts",24,0.4,1.0,40,0.4,2.4);
	h2p_list.push_back(h2p_pmiss_mom1);
	TH2D * h2p_pmiss_mom2 = new TH2D("epp_pmiss_mom2","pmiss_epsilon;pmiss;pRec;Counts",24,0.4,1.0,17,0.35,1.2);
	h2p_list.push_back(h2p_pmiss_mom2);
	TH2D * h2p_pmiss_momrat = new TH2D("epp_pmiss_momrat","pmiss_epsilon;pmiss;pLead/pRec;Counts",24,0.4,1.0,30,1.,4.);
	h2p_list.push_back(h2p_pmiss_momrat);
	TH2D * h2p_pmiss_momdiff = new TH2D("epp_pmiss_momdiff","pmiss_epsilon;pmiss;pLead-pRec;Counts",24,0.4,1.0,40,0.,2.);
	h2p_list.push_back(h2p_pmiss_momdiff);

	TH2D * h2p_pRec_epsilon = new TH2D("epp_pRec_epsilon","pRec_epsilon;pRec;epsilon;Counts",20,0.3,1.2,20,0.2,1.2);
	h2p_list.push_back(h2p_pRec_epsilon);
	TH1D * h2p_pRec_epsilon_mean = new TH1D("epp_pRec_epsilon_mean","pRec_epsilon_mean;pRec;epsilon_mean;Counts",20,0.3,1.2);
	h2p_pRec_epsilon_mean->Sumw2();
	TH1D * h2p_pRec_epsilon_std = new TH1D("epp_pRec_epsilon_std","pRec_epsilon_std;pRec;epsilon_std;Counts",20,0.3,1.2);
	h2p_pRec_epsilon_std->Sumw2();
	TH2D * h2p_pRec_eMiss = new TH2D("epp_pRec_eMiss","pRec_eMiss;pRec;eMiss;Counts",10,0.35,0.9,10,0,0.5);
	h2p_list.push_back(h2p_pRec_eMiss);
	TH1D * h2p_pRec_eMiss_mean = new TH1D("epp_pRec_eMiss_mean","pRec_eMiss_mean;pRec;eMiss_mean;Counts",20,0.3,1.2);
	h2p_pRec_eMiss_mean->Sumw2();
	TH1D * h2p_pRec_eMiss_std = new TH1D("epp_pRec_eMiss_std","pRec_eMiss_std;pRec;eMiss_std;Counts",20,0.3,1.2);
	h2p_pRec_eMiss_std->Sumw2();
	TH1D * h2p_pRecError = new TH1D("epp_pRecError","epp;pRecError;Counts",40,-1,1);
	h2p_pRecError->Sumw2();
	TH2D * pp_to_p_2d = new TH2D("pp_to_p_2d","2d ratio;pmiss [GeV];E1 [GeV];pp/p",28,0.35,1.0,20,0.5,0.9);
	
	TH1D * h1p_Emiss = new TH1D("ep_Emiss","ep;Emiss [GeV];Counts",40,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss);
	TH1D * h1p_Emiss_fine = new TH1D("ep_Emiss_fine","ep;Emiss [GeV];Counts",160,-0.2,0.6);
	h1p_list.push_back(h1p_Emiss_fine);
	TH1D * h1p_e1 = new TH1D("ep_e1","ep;e1 [GeV];Counts",160,0.5,1.0);
	h1p_list.push_back(h1p_e1);
	TH1D * h1p_emiss = new TH1D("ep_emiss","ep;emiss [GeV];Counts",160,0.0,0.5);
	h1p_list.push_back(h1p_emiss);
	TH1D * h2p_Emiss = new TH1D("epp_Emiss","epp;Emiss [GeV];Counts",40,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss);
	TH1D * h2p_Emiss_fine = new TH1D("epp_Emiss_fine","epp;Emiss [GeV];Counts",160,-0.2,0.6);
	h2p_list.push_back(h2p_Emiss_fine);	
	TH1D * h2p_e1 = new TH1D("epp_e1","epp;e1 [GeV];Counts",160,0.5,1.0);
	h2p_list.push_back(h2p_e1);
	TH1D * h2p_emiss = new TH1D("epp_emiss","epp;emiss [GeV];Counts",160,0.0,1.5);
	h2p_list.push_back(h2p_emiss);

	//The first element is pmiss bin, second is xB bin, third is QSq bin
	TH1D * h1p_Emiss_split[4][3][3]; 
	TH1D * h1p_Emiss_fine_split[4][3][3]; 
	TH1D * h1p_e1_split[4][3][3]; 
	TH1D * h1p_emiss_split[4][3][3]; 
	TH1D * h1p_Pmq_split[4][3][3]; 
	TH1D * h1p_Pmzq_split[4][3][3]; 
	TH1D * h1p_PmTq_split[4][3][3];
	TH1D * h1p_dE1_split[4][3][3];
	TH1D * h1p_rE1_split[4][3][3];
	 
	TH1D * h2p_Emiss_split[4][3][3]; 
	TH1D * h2p_Emiss_fine_split[4][3][3];
	TH1D * h2p_e1_split[4][3][3]; 
	TH1D * h2p_emiss_split[4][3][3]; 
	TH1D * h2p_Pmq_split[4][3][3]; 
	TH1D * h2p_Pmzq_split[4][3][3]; 
	TH1D * h2p_PmTq_split[4][3][3]; 
	TH1D * h2p_dE1_split[4][3][3];
	TH1D * h2p_rE1_split[4][3][3];
	
	//dir_sub_bins->cd();
	for (int i=0 ; i<4 ; i++)
	  {
	    for(int j=0; j<3 ; j++){
	  
	      for(int k=0; k<3; k++){
		char temp[100];
		
		sprintf(temp,"ep_Emiss_%d_%d_%d",i,j,k);
		h1p_Emiss_split[i][j][k] = new TH1D(temp,"ep;Emiss [GeV];Counts",40,-0.2,0.6);
		h1p_list.push_back(h1p_Emiss_split[i][j][k]);
		
		sprintf(temp,"ep_Emiss_fine_%d_%d_%d",i,j,k);
		h1p_Emiss_fine_split[i][j][k] = new TH1D(temp,"ep;Emiss [GeV];Counts",160,-0.2,0.6);
		h1p_list.push_back(h1p_Emiss_fine_split[i][j][k]);

		sprintf(temp,"ep_e1_%d_%d_%d",i,j,k);
		h1p_e1_split[i][j][k] = new TH1D(temp,"ep;e1 [GeV];Counts",120,0.4,1.0);
		h1p_list.push_back(h1p_e1_split[i][j][k]);
		
		sprintf(temp,"ep_emiss_%d_%d_%d",i,j,k);
		h1p_emiss_split[i][j][k] = new TH1D(temp,"ep;emiss [GeV];Counts",120,-0.1,0.5);
		h1p_list.push_back(h1p_emiss_split[i][j][k]);
		
		sprintf(temp,"ep_Pmq_%d_%d_%d",i,j,k);
		h1p_Pmq_split[i][j][k] = new TH1D(temp,"ep;Pmq [GeV];Counts",20,100.,180.);
		h1p_list.push_back(h1p_Pmq_split[i][j][k]);
		
		sprintf(temp,"ep_Pmzq_%d_%d_%d",i,j,k);
		h1p_Pmzq_split[i][j][k] = new TH1D(temp,"ep;Pmzq [GeV];Counts",30,-1.0,0.0);
		h1p_list.push_back(h1p_Pmzq_split[i][j][k]);
		
		sprintf(temp,"ep_PmTq_%d_%d_%d",i,j,k);
		h1p_PmTq_split[i][j][k] = new TH1D(temp,"ep;PmTq [GeV];Counts",30,0.0,1.0);
		h1p_list.push_back(h1p_PmTq_split[i][j][k]);

		sprintf(temp,"ep_dE1_%d_%d_%d",i,j,k);
		h1p_dE1_split[i][j][k] = new TH1D(temp,"ep;dE1 [GeV];Counts",40,-0.2,1);
		h1p_list.push_back(h1p_dE1_split[i][j][k]);

		sprintf(temp,"ep_rE1_%d_%d_%d",i,j,k);
		h1p_rE1_split[i][j][k] = new TH1D(temp,"ep;rE1;Counts",40,0,1);
		h1p_list.push_back(h1p_rE1_split[i][j][k]);
		
		sprintf(temp,"epp_Emiss_%d_%d_%d",i,j,k);
		h2p_Emiss_split[i][j][k] = new TH1D(temp,"epp;Emiss [GeV];Counts",20,-0.2,0.6);
		h2p_list.push_back(h2p_Emiss_split[i][j][k]);
		
		sprintf(temp,"epp_Emiss_fine_%d_%d_%d",i,j,k);
		h2p_Emiss_fine_split[i][j][k] = new TH1D(temp,"epp;Emiss [GeV];Counts",80,-0.2,0.6);
		h2p_list.push_back(h2p_Emiss_fine_split[i][j][k]);

		sprintf(temp,"epp_e1_%d_%d_%d",i,j,k);
		h2p_e1_split[i][j][k] = new TH1D(temp,"epp;e1 [GeV];Counts",60,0.4,1.0);
		h2p_list.push_back(h2p_e1_split[i][j][k]);

		sprintf(temp,"epp_emiss_%d_%d_%d",i,j,k);
		h2p_emiss_split[i][j][k] = new TH1D(temp,"epp;emiss [GeV];Counts",60,-0.1,0.5);
		h2p_list.push_back(h2p_emiss_split[i][j][k]);

		sprintf(temp,"epp_Pmq_%d_%d_%d",i,j,k);
		h2p_Pmq_split[i][j][k] = new TH1D(temp,"epp;Pmq [GeV];Counts",20,100.,180.);
		h2p_list.push_back(h2p_Pmq_split[i][j][k]);
		
		sprintf(temp,"epp_Pmzq_%d_%d_%d",i,j,k);
		h2p_Pmzq_split[i][j][k] = new TH1D(temp,"epp;Pmzq [GeV];Counts",30,-1.0,0.0);
		h2p_list.push_back(h2p_Pmzq_split[i][j][k]);

		sprintf(temp,"epp_PmTq_%d_%d_%d",i,j,k);
		h2p_PmTq_split[i][j][k] = new TH1D(temp,"epp;PmTq [GeV];Counts",30,0.0,1.0);
		h2p_list.push_back(h2p_PmTq_split[i][j][k]);

		sprintf(temp,"epp_dE1_%d_%d_%d",i,j,k);
		h2p_dE1_split[i][j][k] = new TH1D(temp,"epp;dE1 [GeV];Counts",40,-0.2,1);
		h2p_list.push_back(h2p_dE1_split[i][j][k]);

		sprintf(temp,"epp_rE1_%d_%d_%d",i,j,k);
		h2p_rE1_split[i][j][k] = new TH1D(temp,"epp;rE1;Counts",40,0,1);
		h2p_list.push_back(h2p_rE1_split[i][j][k]); 

	      }
	    }
	  }

	dir_by_sec->cd();
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

	fo->cd();

	// Now that all histograms have been defined, set them to Sumw2
	for (int i=0 ; i<h1p_list.size() ; i++)
	  h1p_list[i]->Sumw2();
	for (int i=0 ; i<h2p_list.size() ; i++)
	  h2p_list[i]->Sumw2();

	// For data and bin-centering
	TH1D * h1p_Pm_30bin =  new TH1D("ep_Pm_30bin" ,"ep;pMiss [GeV];Counts",30,0.4,1.0);
	h1p_Pm_30bin->Sumw2();
	TH1D * h2p_Pm_30bin =  new TH1D("epp_Pm_30bin" ,"epp;pMiss [GeV];Counts",30,0.4,1.0);
	h2p_Pm_30bin->Sumw2();
	TH1D * h1p_Pm_30bin_bins =  new TH1D("ep_Pm_30bin_bins" ,"ep;pMiss [GeV];Sum pMiss [GeV]",30,0.4,1.0);
	TH1D * h2p_Pm_30bin_bins =  new TH1D("epp_Pm_30bin_bins" ,"epp;pMiss [GeV];Sum pMiss [GeV]",30,0.4,1.0);
	TH1D * h1p_Pm_coarse_bins =  new TH1D("ep_Pm_coarse_bins" ,"ep;pMiss [GeV];Sum pMiss [GeV]",9,coarse_bin_edges_new);
	TH1D * h2p_Pm_coarse_bins =  new TH1D("epp_Pm_coarse_bins" ,"epp;pMiss [GeV];Sum pMiss [GeV]",9,coarse_bin_edges_new);

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
	bool resetto1 = false;
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
	else resetto1 = true;
	
	for (int event =0 ; event < t1p->GetEntries() ; event++)
	{
		t1p->GetEvent(event);

		if(resetto1) weight = 1;

		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);

		if (doOtherCut) {
		  // Do necessary cuts
		  if (fabs(Rp[0][2]+22.25)>2.25)
		    continue;
		  if (Pp_size[0]>2.4)
		    continue;
		  if (Pmiss_size[0]<pmiss_cut)
		    continue;
		}
		
		if (doCut){
		  // Apply fiducial cuts
		  if (!accept_electron(ve))
		    continue;
		  if (doGaps)
		    {
		      if (!accept_proton(vp))
			continue;
		    }
		  else
		    {
		      if (!accept_proton_simple(vp))
			continue;
		    }

		  // Apply an additional fiducial cut that the map acc must be > acc_thresh
		  if ( proton_map.accept(vp) < acc_thresh)
		    continue;
		}

		// Sector-specific theta1 cuts
		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vp.Theta() * 180./M_PI;
		
                // A few more vectors                                                            
                TVector3 vq(q[0],q[1],q[2]);
                TVector3 vqUnit = vq.Unit();
		TVector3 vm = vp - vq;

		double omega = Q2/(2.*mN*Xb);

 		//Do spectral function weight
		if (doSWeight){
		  double factorS = omega * myCS.sigma_eN(Ebeam,ve,vp,true) / (2 * Ebeam * ve.Mag() * Xb);
		  if (factorS>sCutOff){
		    weight = weight / factorS;  
		  }
		  else weight=0;
		}
				
		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pm_30bin ->Fill(Pmiss_size[0],weight);
		h1p_Pm_30bin_bins ->Fill(Pmiss_size[0],weight*Pmiss_size[0]);
		h1p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h1p_Pm_coarse_bins->Fill(Pmiss_size[0],weight*Pmiss_size[0]);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		// Kinematic variables we need
		double Ep = sqrt(Pp_size[0]*Pp_size[0] + mN*mN);
		double Emiss = -m_12C + mN + sqrt( sq(omega + m_12C - Ep) - (Pmiss_size[0]*Pmiss_size[0]));
		double epsilon = Ep - omega;
		double dE1 = sqrt(vm.Mag2()+sq(mN)) - epsilon;
		double rE1 = epsilon/sqrt(vm.Mag2()+sq(mN));

		//Let's calculate light cone variables                                           
                double alphaq= (omega - vq.Mag()) / mN;
                double alphaLead = (Ep - vp.Dot(vqUnit)) / mN;
                double alphaM = alphaLead - alphaq;

                h1p_alphaq->Fill(alphaq,weight);
                h1p_alphaLead->Fill(alphaLead,weight);
                h1p_alphaM->Fill(alphaM,weight);
		h1p_dE1->Fill(dE1,weight);
		h1p_rE1->Fill(rE1,weight);
		
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

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		
		// Let's figure out missing energy! 
		h1p_Emiss->Fill(Emiss,weight);
		h1p_Emiss_fine->Fill(Emiss,weight);
		h1p_e1->Fill(epsilon,weight);
		h1p_emiss->Fill(mN-epsilon,weight);
		h1p_pmiss_Emiss->Fill(Pmiss_size[0],Emiss,weight);
		h1p_Emiss_by_sector->Fill(sec_e,Emiss);

		//Pmiss splits
		int Pmiss_region_p;
	       	if (Pmiss_size[0] < pmiss_lo) Pmiss_region_p = 0;
		else if (Pmiss_size[0] < pmiss_md) Pmiss_region_p = 1;
		else if (Pmiss_size[0] < pmiss_hi) Pmiss_region_p = 2;
		else Pmiss_region_p = 3;
		int xB_region_p;
		if (Xb < 1.4) xB_region_p = 1;
		else xB_region_p = 2;
		int QSq_region_p;
		if (Q2 < 2) QSq_region_p = 1;
		else  QSq_region_p = 2;

		for(int j=0; j<=xB_region_p; j=j+xB_region_p){
		  for(int k=0; k<=QSq_region_p; k=k+QSq_region_p){
		    h1p_Emiss_split[Pmiss_region_p][j][k]->Fill(Emiss,weight);
		    h1p_Emiss_fine_split[Pmiss_region_p][j][k]->Fill(Emiss,weight);
		    h1p_e1_split[Pmiss_region_p][j][k]->Fill(epsilon,weight);
		    h1p_emiss_split[Pmiss_region_p][j][k]->Fill(mN-epsilon,weight);
		    h1p_Pmq_split[Pmiss_region_p][j][k]->Fill(Pmiss_q_angle[0],weight);
		    h1p_Pmzq_split[Pmiss_region_p][j][k]->Fill(vm.Dot(vqUnit),weight);
		    h1p_PmTq_split[Pmiss_region_p][j][k]->Fill(vm.Perp(vqUnit),weight);
		    h1p_dE1_split[Pmiss_region_p][j][k]->Fill(dE1,weight);
		    h1p_rE1_split[Pmiss_region_p][j][k]->Fill(rE1,weight);
		  }
		}
		
		h1p_pmiss_E1->Fill(Pmiss_size[0],epsilon,weight);
		h1p_pmiss_epsilon->Fill(Pmiss_size[0],epsilon,weight);

		h1p_pmiss_QSq->Fill(Pmiss_size[0],Q2,weight);
		h1p_pmiss_mom1->Fill(Pmiss_size[0],Pp_size[0],weight);
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
	resetto1 = false;
	weight_branch = t2p->GetBranch("weight");
	if (weight_branch)
	{
		t2p->SetBranchAddress("weight",&weight);
	}
	else resetto1 = true;
	for (int event =0 ; event < t2p->GetEntries() ; event++)
	{
		t2p->GetEvent(event);

		if(resetto1) weight = 1;

		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vlead(Pp[0][0],Pp[0][1],Pp[0][2]);
	
		if(doOtherCut){
		  // Do necessary cuts
		  if (fabs(Rp[0][2]+22.25)>2.25)
		    continue;
		  if (Pp_size[0]>2.4)
		    continue;
		  if (Pmiss_size[0]<pmiss_cut)
		    continue;
		}
		
		if(doCut){
		  // Apply fiducial cuts
		  if (!accept_electron(ve))
		    continue;
		  if (doGaps)
		    {
		      if (!accept_proton(vlead))
			continue;
		    }
		  else
		    {
		      if (!accept_proton_simple(vlead))
			continue;
		    }

		  // Apply an additional fiducial cut that the map acc must be > acc_thresh
		  if ( proton_map.accept(vlead) < acc_thresh)
		    continue;
		}

		// Sector-specific theta1 cuts
		double phi1_deg = vlead.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vlead.Theta() * 180./M_PI;

                // A few more vectors                                                       
                TVector3 vq(q[0],q[1],q[2]);
                TVector3 vqUnit = vq.Unit();
                TVector3 vmiss = vlead - vq;
		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);      

		double omega = Q2/(2.*mN*Xb);

		//Do spectral function weight
		if (doSWeight){
		  double factorS = omega * myCS.sigma_eN(Ebeam,ve,vlead,true) / (2 * Ebeam * ve.Mag() * Xb);
		  if (factorS>sCutOff){
		    weight = weight / factorS;
		  }
		  else weight=0;
		}
		
		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_Pm_30bin ->Fill(Pmiss_size[0],weight);
		h1p_Pm_30bin_bins ->Fill(Pmiss_size[0],weight*Pmiss_size[0]);
		h1p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h1p_Pm_coarse_bins->Fill(Pmiss_size[0],weight*Pmiss_size[0]);
		h1p_Pmq->Fill(Pmiss_q_angle[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		// Kinematic variables we need
		double Elead = sqrt(Pp_size[0]*Pp_size[0] + mN*mN);
		double Emiss = -m_12C + mN + sqrt( sq(omega + m_12C - Elead) - (Pmiss_size[0]*Pmiss_size[0]));
		double epsilon = Elead - omega;
		double dE1 = sqrt(vmiss.Mag2()+sq(mN)) - epsilon;
		double rE1 = epsilon/sqrt(vmiss.Mag2()+sq(mN));

		// Calculate the expected recoil momentum
		double TCM = (3/20) * sq(sigmaCM) / mN;
		double p2Calc = sqrt( sq( (m_12C-(m_10B+Estar)) - epsilon - TCM ) - sq(mN));
		double p2CalcError = (p2Calc - vrec.Mag()) / vrec.Mag();

		h2p_pRecError->Fill(p2CalcError,weight);
		
                // Let's calculate light cone variables                                                      
                double alphaq= (omega - vq.Mag()) / mN;
                double alphaLead = (Elead - vlead.Dot(vqUnit)) / mN;
                double alphaM = alphaLead - alphaq;

                h2p_alphaq->Fill(alphaq,weight);
                h2p_alphaLead->Fill(alphaLead,weight);
                h2p_alphaM->Fill(alphaM,weight);
		h2p_dE1->Fill(dE1,weight);
		h2p_rE1->Fill(rE1,weight);
		
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

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		
		h1p_Emiss->Fill(Emiss,weight);
		h1p_Emiss_fine->Fill(Emiss,weight);
		h1p_e1->Fill(epsilon,weight);
		h1p_emiss->Fill(mN-epsilon,weight);
		h1p_pmiss_Emiss->Fill(Pmiss_size[0],Emiss,weight);
		h1p_pmiss_E1->Fill(Pmiss_size[0],epsilon,weight);
		h1p_Emiss_by_sector->Fill(sec_e,Emiss);


		//Look at kinematic regions
		int Pmiss_region_pp;
	       	if (Pmiss_size[0] < pmiss_lo) Pmiss_region_pp = 0;
		else if (Pmiss_size[0] < pmiss_md) Pmiss_region_pp = 1;
		else if (Pmiss_size[0] < pmiss_hi) Pmiss_region_pp = 2;
		else Pmiss_region_pp = 3;
		int xB_region_pp;
		if (Xb < 1.4) xB_region_pp = 1;
		else xB_region_pp = 2;
		int QSq_region_pp;
		if (Q2 < 2) QSq_region_pp = 1;
		else  QSq_region_pp = 2;
		

		for(int j=0; j<=xB_region_pp; j=j+xB_region_pp){
		  for(int k=0; k<=QSq_region_pp; k=k+QSq_region_pp){
		    h1p_Emiss_split[Pmiss_region_pp][j][k]->Fill(Emiss,weight);
		    h1p_Emiss_fine_split[Pmiss_region_pp][j][k]->Fill(Emiss,weight);
		    h1p_e1_split[Pmiss_region_pp][j][k]->Fill(epsilon,weight);
		    h1p_emiss_split[Pmiss_region_pp][j][k]->Fill(mN-epsilon,weight);
		    h1p_Pmq_split[Pmiss_region_pp][j][k]->Fill(Pmiss_q_angle[0],weight);
		    h1p_Pmzq_split[Pmiss_region_pp][j][k]->Fill(vmiss.Dot(vqUnit),weight);
		    h1p_PmTq_split[Pmiss_region_pp][j][k]->Fill(vmiss.Perp(vqUnit),weight);
		    h1p_dE1_split[Pmiss_region_pp][j][k]->Fill(dE1,weight);
		    h1p_rE1_split[Pmiss_region_pp][j][k]->Fill(rE1,weight);
		  }
		}
		h1p_pmiss_epsilon->Fill(Pmiss_size[0],epsilon,weight);
		
		h1p_pmiss_QSq->Fill(Pmiss_size[0],Q2,weight);
		h1p_pmiss_mom1->Fill(Pmiss_size[0],Pp_size[0],weight);
		

		if(doOtherCut){
		  // Make a check on the recoils
		  if (fabs(Rp[1][2]+22.25)>2.25)
		    continue;
		  if (Pp_size[1] < 0.35)
		    continue;
		}
		
		if(doCut){

		  if (doGaps)
		    {
		      if (!accept_proton(vrec))
			continue;
		    }
		  else
		    {
		      if (!accept_proton_simple(vrec))
			continue;
		    }


		  // Apply an additional fiducial cut that the map acc must be > acc_thresh
		  if ( proton_map.accept(vrec) < acc_thresh)
		    continue;
		}

                double Erec = sqrt(vrec.Mag2() + mN*mN);
		TVector3 vcm = vmiss + vrec;

                // Find final light cone variables                                                          
                double alphaRec = (Erec - vrec.Dot(vqUnit))/mN;
                double alphaD = alphaM + alphaRec;

		h2p_QSq->Fill(Q2,weight);
		h2p_xB ->Fill(Xb,weight);
		h2p_Pm ->Fill(Pmiss_size[0],weight);
		h2p_Pm_30bin ->Fill(Pmiss_size[0],weight);
		h2p_Pm_30bin_bins ->Fill(Pmiss_size[0],Pmiss_size[0]*weight);
		h2p_Pm_clas ->Fill(Pmiss_size[0],weight);
		h2p_Pm_coarse->Fill(Pmiss_size[0],weight);
		h2p_Pm_coarse_bins->Fill(Pmiss_size[0],weight*Pmiss_size[0]);
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
		h2p_e1->Fill(epsilon,weight);
		h2p_emiss->Fill(mN-epsilon,weight);
		h2p_pmiss_E1->Fill(Pmiss_size[0],epsilon,weight);

		h2p_pmiss_epsilon->Fill(Pmiss_size[0],epsilon,weight);
		h2p_pRec_epsilon->Fill(Pp_size[1],epsilon,weight);
	        h2p_pRec_eMiss->Fill(Pp_size[1],Emiss,weight);

		//Do Pmiss splits
		
		for(int j=0; j<=xB_region_pp; j=j+xB_region_pp){
		  for(int k=0; k<=QSq_region_pp; k=k+QSq_region_pp){
		    h2p_Emiss_split[Pmiss_region_pp][j][k]->Fill(Emiss,weight);
		    h2p_Emiss_fine_split[Pmiss_region_pp][j][k]->Fill(Emiss,weight);
		    h2p_e1_split[Pmiss_region_pp][j][k]->Fill(epsilon,weight);
		    h2p_emiss_split[Pmiss_region_pp][j][k]->Fill(mN-epsilon,weight);
		    h2p_Pmq_split[Pmiss_region_pp][j][k]->Fill(Pmiss_q_angle[0],weight);
		    h2p_Pmzq_split[Pmiss_region_pp][j][k]->Fill(vmiss.Dot(vqUnit),weight);
		    h2p_PmTq_split[Pmiss_region_pp][j][k]->Fill(vmiss.Perp(vqUnit),weight);
		    h2p_dE1_split[Pmiss_region_pp][j][k]->Fill(dE1,weight);
		    h2p_rE1_split[Pmiss_region_pp][j][k]->Fill(rE1,weight);
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

		h2p_pmiss_QSq->Fill(Pmiss_size[0],Q2,weight);
		h2p_pmiss_mom1->Fill(Pmiss_size[0],Pp_size[0],weight);
		h2p_pmiss_mom2->Fill(Pmiss_size[0],Pp_size[1],weight);
		h2p_pmiss_momrat->Fill(Pmiss_size[0],Pp_size[0]/Pp_size[1],weight);

		TVector3 vdiff = vlead - vrec;

		h2p_momdiff->Fill(vdiff.Mag(),weight);
		h2p_pmiss_momdiff->Fill(Pmiss_size[0],vdiff.Mag(),weight);

		// Histogram for the "apparent E*"
		double apparent_Estar = sqrt(sq(sqrt(sq(m_10B) + vcm.Mag2()) + Erec) -vlead.Mag2()) - m_11B;
		h2p_pmiss_appEstar->Fill(Pmiss_size[0],apparent_Estar,weight);
	}

	for( int i = 0; i < h2p_pRec_epsilon->GetNbinsX(); i++){
	  TH1D * h2p_pRec_epsilon_Proj = h2p_pRec_epsilon->ProjectionY("pRec_epsilon_Proj",i,i+1);
	  h2p_pRec_epsilon_mean->Fill(h2p_pRec_epsilon->GetXaxis()->GetBinCenter(i+1),h2p_pRec_epsilon_Proj->GetMean());
	  h2p_pRec_epsilon_std->Fill(h2p_pRec_epsilon->GetXaxis()->GetBinCenter(i+1),h2p_pRec_epsilon_Proj->GetStdDev());
	}
	for( int i = 0; i < h2p_pRec_eMiss->GetNbinsX(); i++){
	  TH1D * h2p_pRec_eMiss_Proj = h2p_pRec_eMiss->ProjectionY("pRec_eMiss_Proj",i,i+1);
	  h2p_pRec_eMiss_mean->Fill(h2p_pRec_eMiss->GetXaxis()->GetBinCenter(i+1),h2p_pRec_eMiss_Proj->GetMean());
	  h2p_pRec_eMiss_std->Fill(h2p_pRec_eMiss->GetXaxis()->GetBinCenter(i+1),h2p_pRec_eMiss_Proj->GetStdDev());
	}	

	f1p->Close();
	f2p->Close();

	// Do the bin centering
	TGraphAsymmErrors * g1p_Pm = new TGraphAsymmErrors(30);
	g1p_Pm->SetName("ep_Pm_graph");
	g1p_Pm->SetTitle("ep;p_miss [GeV];Counts");
	for (int i=1 ; i<= 30 ; i++)
	  {
	    // e'p
	    double sumN = h1p_Pm_30bin->GetBinContent(i);
	    double sumPm = h1p_Pm_30bin_bins->GetBinContent(i);
	    double avgPm = h1p_Pm_30bin->GetBinLowEdge(i) + 0.5 *h1p_Pm_30bin->GetBinWidth(i);
	    if (sumN > 0) avgPm = (sumPm/sumN);
	    h1p_Pm_30bin_bins->SetBinContent(i,avgPm);

	    // Fill out the TGraph
	    g1p_Pm->SetPoint(i-1,avgPm,sumN);
	    g1p_Pm->SetPointError(i-1,avgPm - h1p_Pm_30bin->GetBinLowEdge(i),
					h1p_Pm_30bin->GetBinLowEdge(i) + h1p_Pm_30bin->GetBinWidth(i) - avgPm,
					h1p_Pm_30bin->GetBinError(i), h1p_Pm_30bin->GetBinError(i));
	  }
	g1p_Pm->Write();
	
	TGraphAsymmErrors * g2p_Pm = new TGraphAsymmErrors(24);
	g2p_Pm->SetName("epp_Pm_graph");
	g2p_Pm->SetTitle("epp;p_miss [GeV];Counts");
	// The binning is non-uniform, so we have to be careful
	const int startBins[24]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,25,27};
	const int stopBins[24] ={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,30};
	for (int i=0 ; i< 24 ; i++)
	  {
	    // e'pp
	    double sumN = h2p_Pm_30bin->Integral(startBins[i],stopBins[i]);
	    double sumPm = h2p_Pm_30bin_bins->Integral(startBins[i],stopBins[i]);
	    double avgPm = 0.5*(h2p_Pm_30bin->GetBinLowEdge(startBins[i]) + h2p_Pm_30bin->GetXaxis()->GetBinUpEdge(stopBins[i]));
	    if (sumN > 0) avgPm = (sumPm / sumN);

	    // Fill out the TGraph
	    //	    cout << startBins[i] << " " << stopBins[i] << "       " 
	    //	 << avgPm << " - " << avgPm - h2p_Pm_30bin->GetBinLowEdge(startBins[i]) << " + "
	    // << h2p_Pm_30bin->GetXaxis()->GetBinUpEdge(stopBins[i]) -avgPm << "         "
	    // << sumN << " - " << sqrt(sumN) << "  + " << sumN << "\n";
	    
	    double binScale = (0.02)/(h2p_Pm_30bin->GetXaxis()->GetBinUpEdge(stopBins[i]) - h2p_Pm_30bin->GetBinLowEdge(startBins[i]));
	    g2p_Pm->SetPoint(i,avgPm,sumN * binScale);
	    g2p_Pm->SetPointError(i,avgPm - h2p_Pm_30bin->GetBinLowEdge(startBins[i]),
				  h2p_Pm_30bin->GetXaxis()->GetBinUpEdge(stopBins[i]) -avgPm,
				  binScale * sqrt(sumN),binScale * sqrt(sumN));
	  }
	g2p_Pm->Write();

	cerr << "The ep and epp integrals are: " << h1p_Pm->Integral() << " "  << h2p_Pm->Integral() << "\n";
	cerr << "Broken down by pmiss range...\n\n";
	for (int j=1 ; j<=30 ; j+=1)
	  {
	    cerr << h1p_Pm->GetBinCenter(j) << " " 
		 << h1p_Pm->GetBinContent(j)
		 << " " 
		 << h2p_Pm->GetBinContent(j)
		 << "\n"; 
	  }
	cerr << "\n\n";

	// pp-to-p, including bin centering
	pp_to_p->BayesDivide(h2p_Pm,h1p_Pm);
	pp_to_p_coarse->BayesDivide(h2p_Pm_coarse,h1p_Pm_coarse);
	for (int j=1 ; j<=9 ; j++)
	  {
	    double sumPm = h1p_Pm_coarse_bins->GetBinContent(j);
	    double sumN = h1p_Pm_coarse->GetBinContent(j);
	    double avgPm = sumPm / sumN;
	    double pp2p = pp_to_p_coarse->GetY()[j-1];
	    pp_to_p_coarse->SetPoint(j-1,avgPm,pp2p);
	    pp_to_p_coarse->SetPointEXlow(j-1,avgPm - h1p_Pm_coarse->GetBinLowEdge(j));
	    pp_to_p_coarse->SetPointEXhigh(j-1,h1p_Pm_coarse->GetXaxis()->GetBinUpEdge(j) - avgPm);

	    cerr << h1p_Pm_coarse->GetBinContent(j)
		 << " & " 
		 << h2p_Pm_coarse->GetBinContent(j)
		 << " & "
		 << pp_to_p_coarse->GetY()[j-1]
		 << " & \n";
	  }
	cerr << "\n\n";

	// 2d pp-to-p
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

	const double data_ep = 5604.;
	const double data_ep_cor = 6077.;
	const double data_epp = 364.;
	const double pnorm = data_ep/h1p_Pm->Integral();
	const double ppnorm = data_epp/h2p_Pm->Integral();

	h2p_pRecError->Scale(data_epp/h2p_Pm->Integral());
	
	// Including a factor if we want to rescale data to match epp luminosity
	const double renorm = data_epp/h2p_Pm->Integral()/(data_ep/h1p_Pm->Integral());
	TVectorT<double> renorm_factor(1);
	renorm_factor[0] = renorm;
	renorm_factor.Write("factor");

	h2p_Pm_clas->Scale(data_epp/h2p_Pm->Integral());

	// scale all the histograms
	for (int i=0 ; i<h1p_list.size() ; i++)
	  h1p_list[i]->Scale(pnorm);
	for (int i=0 ; i<h2p_list.size() ; i++)
	  h2p_list[i]->Scale(ppnorm);

	// Write all the histograms and close the file.
	fo->Write();
	fo->Close();

	return 0;
}

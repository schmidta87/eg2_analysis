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

const double pmiss_lo=0.5;
const double pmiss_md=0.6;
const double pmiss_hi=0.7;

const double pmiss_cut=0.4;

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

	// Create histograms
	TH1D * h1p_QSq = new TH1D("ep_QSq","ep;QSq [GeV^2];Counts",40,1.,5.);
	h1p_list.push_back(h1p_QSq);
	TH1D * h1p_xB =  new TH1D("ep_xB" ,"ep;xB;Counts",26,1.2,2.5);
	h1p_list.push_back(h1p_xB );
	TH1D * h1p_thetam = new TH1D("ep_thetam","ep;Theta_m [deg];Counts",40,10.,170.);
	h1p_list.push_back(h1p_thetam);
	TH1D * h1p_cthetam = new TH1D("ep_cthetam","ep;cos Theta_m;Counts",40,-1,1);
	h1p_list.push_back(h1p_cthetam);
	TH1D * h1p_Pmz =  new TH1D("ep_Pmz" ,"ep;pMiss_z [GeV];Counts",30,0.,0.4);
	h1p_list.push_back(h1p_Pmz);
	TH1D * h1p_PmT =  new TH1D("ep_PmT" ,"ep;pMiss_T [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmT);
	TH1D * h1p_Pmzq =  new TH1D("ep_Pmzq" ,"ep;pMiss_zq [GeV];Counts",30,-1.0,0.0);
	h1p_list.push_back(h1p_Pmzq);
	TH1D * h1p_PmTq =  new TH1D("ep_PmTq" ,"ep;pMiss_Tq [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h1p_PmTq);
	TH1D * h2p_Pmzrel =  new TH1D("epp_Pmzrel" ,"ep;pMiss_zrel [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h2p_Pmzrel);
	TH1D * h2p_PmTrel =  new TH1D("epp_PmTrel" ,"ep;pMiss_Trel [GeV];Counts",30,0.0,0.5);
	h1p_list.push_back(h2p_PmTrel);
	TH1D * h2p_Pmzrec =  new TH1D("epp_Pmzrec" ,"ep;pMiss_zrec [GeV];Counts",30,-1.0,1.0);
	h1p_list.push_back(h2p_Pmzrec);
	TH1D * h2p_PmTrec =  new TH1D("epp_PmTrec" ,"ep;pMiss_Trec [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h2p_PmTrec);
	TH1D * h2p_PmzCM =  new TH1D("epp_PmzCM" ,"ep;pMiss_zCM [GeV];Counts",30,-1.0,1.0);
	h1p_list.push_back(h2p_PmzCM);
	TH1D * h2p_PmTCM =  new TH1D("epp_PmTCM" ,"ep;pMiss_TCM [GeV];Counts",30,0.0,1.0);
	h1p_list.push_back(h2p_PmTCM);
	TH1D * h1p_Pm =  new TH1D("ep_Pm" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm);
	TH1D * h2p_Pm =  new TH1D("epp_Pm" ,"epp;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h2p_Pm);
	TH1D * h1p_Pm_10 =  new TH1D("ep_Pm_10" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_10);
	TH1D * h1p_Pm_15 =  new TH1D("ep_Pm_15" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_15);
	TH1D * h1p_Pm_17 =  new TH1D("ep_Pm_17" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_17);
	TH1D * h1p_Pm_20 =  new TH1D("ep_Pm_20" ,"ep;pMiss [GeV];Counts",28,0.3,1.0);
	h1p_list.push_back(h1p_Pm_20);
	TH1D * h1p_cPeq = new TH1D("ep_cPeq","ep;cos(Theta_Peq);Counts",35,0.1,0.45);
	h1p_list.push_back(h1p_cPeq);
	TH1D * h1p_cPmq = new TH1D("ep_cPmq","ep;cos(Theta_Pmq);Counts",40,-1.,0.);
	h1p_list.push_back(h1p_cPmq);
	TH1D * h1p_cP1q = new TH1D("ep_cP1q","ep;cos(Theta_P1q);Counts",30,0.85,1.);
	h1p_list.push_back(h1p_cP1q);
	TH1D * h1p_cPm1 = new TH1D("ep_cPm1","ep;cos(Theta_Pm1);Counts",40,-1.,0.);
	h1p_list.push_back(h1p_cPm1);
	TH1D * h2p_cPm2 = new TH1D("epp_cPm2","ep;cos(Theta_Pm2);Counts",40,-1.,0.);
	h1p_list.push_back(h2p_cPm2);
	TH1D * h1p_theta1 = new TH1D("ep_theta1","ep;Theta_1 [deg];Counts",40,10.,90.);
	h1p_list.push_back(h1p_theta1);
	TH1D * h2p_theta1 = new TH1D("epp_theta1","epp;Theta_1 [deg];Counts",40,10.,90.);
	h1p_list.push_back(h2p_theta1);
	TH1D * h1p_thetae = new TH1D("ep_thetae","ep;Theta_e [deg];Counts",60,10.,40.);
	h1p_list.push_back(h1p_thetae);
	TH1D * h1p_mome = new TH1D("ep_mome","ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
	h1p_list.push_back(h1p_mome);
	TH1D * h1p_momez = new TH1D("ep_momez","ep;Mom_ez [GeV/c];Counts",60,2.7,4.7);
	h1p_list.push_back(h1p_momez);
	TH1D * h1p_momeT = new TH1D("ep_momeT","ep;Mom_eT [GeV/c];Counts",60,1.0,2.0);
	h1p_list.push_back(h1p_momeT);
	TH1D * h1p_momezq = new TH1D("ep_momezq","ep;Mom_ezq [GeV/c];Counts",60,0.0,2.0);
	h1p_list.push_back(h1p_momezq);
	TH1D * h1p_momeTq = new TH1D("ep_momeTq","ep;Mom_eTq [GeV/c];Counts",60,2.5,5.0);
	h1p_list.push_back(h1p_momeTq);
	TH1D * h1p_mom1 = new TH1D("ep_mom1","ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
	h1p_list.push_back(h1p_mom1);
	TH1D * h1p_mom1z = new TH1D("ep_mom1z","ep;Mom_1z [GeV/c];Counts",60,0.0,2.0);
	h1p_list.push_back(h1p_mom1z);
	TH1D * h1p_mom1T = new TH1D("ep_mom1T","ep;Mom_1T [GeV/c];Counts",60,0.0,2.0);
	h1p_list.push_back(h1p_mom1T);
	TH1D * h1p_mom1zq = new TH1D("ep_mom1zq","ep;Mom_1zq [GeV/c];Counts",60,0.5,3.0);
	h1p_list.push_back(h1p_mom1zq);
	TH1D * h1p_mom1Tq = new TH1D("ep_mom1Tq","ep;Mom_1Tq [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h1p_mom1Tq);
	TH1D * h1p_phi1 = new TH1D("ep_phi1","ep;Phi_1 [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phi1);
	TH1D * h1p_phie = new TH1D("ep_phie","ep;Phi_e [deg];Counts",60,-30.,330.);
	h1p_list.push_back(h1p_phie);
	TH1D * h1p_q = new TH1D("ep_q","ep;q [GeV/c];Counts",40,1.0,3.0);
	h1p_list.push_back(h1p_q);
	TH1D * h1p_qz = new TH1D("ep_qz","ep;q_z [GeV/c];Counts",50,0.0,2.0);
	h1p_list.push_back(h1p_qz);
	TH1D * h1p_qT = new TH1D("ep_qT","ep;q_T [GeV/c];Counts",40,1.0,2.0);
	h1p_list.push_back(h1p_qT);
	TH1D * h1p_thetaq = new TH1D("ep_thetaq","ep;Theta_q [deg];Counts",50,30.,80.);
	h1p_list.push_back(h1p_thetaq);	
	TH1D * h2p_mom2 = new TH1D("epp_mom2","epp;Recoil Mom [GeV/c];Counts",17,0.35,1.2);
	h1p_list.push_back(h2p_mom2);
	TH1D * h2p_mom2z = new TH1D("epp_mom2z","ep;Mom_1z [GeV/c];Counts",60,0-1.2,1.2);
	h1p_list.push_back(h2p_mom2z);
	TH1D * h2p_mom2T = new TH1D("epp_mom2T","ep;Mom_1T [GeV/c];Counts",60,00,1.2);
	h1p_list.push_back(h2p_mom2T);
	TH1D * h2p_mom2zq = new TH1D("epp_mom2zq","ep;Mom_1zq [GeV/c];Counts",60,-1.2,1.2);
	h1p_list.push_back(h2p_mom2zq);
	TH1D * h2p_mom2Tq = new TH1D("epp_mom2Tq","ep;Mom_1Tq [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_mom2Tq);
	TH1D * h2p_pCM = new TH1D("epp_pCM","ep;p_CM [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_pCM);
	TH1D * h2p_pCMz = new TH1D("epp_pCMz","ep;p_CMz [GeV/c];Counts",60,-1.0,1.0);
	h1p_list.push_back(h2p_pCMz);
	TH1D * h2p_pCMT = new TH1D("epp_pCMT","ep;p_CMT [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_pCMT);
	TH1D * h2p_pCMzPm = new TH1D("epp_pCMzPm","ep;p_CMzPm [GeV/c];Counts",60,-1.0,1.0);
	h1p_list.push_back(h2p_pCMzPm);
	TH1D * h2p_pCMTPm = new TH1D("epp_pCMTPm","ep;p_CMTPm [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_pCMTPm);
	TH1D * h2p_pCMzq = new TH1D("epp_pCMzq","ep;p_CMzq [GeV/c];Counts",60,-1.0,1.0);
	h1p_list.push_back(h2p_pCMzq);
	TH1D * h2p_pCMTq = new TH1D("epp_pCMTq","ep;p_CMTq [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_pCMTq);
	TH1D * h2p_pCMzPr = new TH1D("epp_pCMzPr","ep;p_CMzPr [GeV/c];Counts",60,-1.0,1.0);
	h1p_list.push_back(h2p_pCMzPr);
	TH1D * h2p_pCMTPr = new TH1D("epp_pCMTPr","ep;p_CMTPr [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_pCMTPr);
	TH1D * h2p_prelz = new TH1D("epp_prelz","ep;p_relz [GeV/c];Counts",60,-1.0,1.0);
	h1p_list.push_back(h2p_prelz);
	TH1D * h2p_prelT = new TH1D("epp_prelT","ep;p_relT [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_prelT);
	TH1D * h2p_prelzPm = new TH1D("epp_prelzPm","ep;p_relzPm [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_prelzPm);
	TH1D * h2p_prelTPm = new TH1D("epp_prelTPm","ep;p_relTPm [GeV/c];Counts",60,0.0,0.7);
	h1p_list.push_back(h2p_prelTPm);
	TH1D * h2p_prelzq = new TH1D("epp_prelzq","ep;p_relzq [GeV/c];Counts",60,-1.0,0.0);
	h1p_list.push_back(h2p_prelzq);
	TH1D * h2p_prelTq = new TH1D("epp_prelTq","ep;p_relTq [GeV/c];Counts",60,0.0,0.7);
	h1p_list.push_back(h2p_prelTq);
	TH1D * h2p_prelzPr = new TH1D("epp_prelzPr","ep;p_relzPr [GeV/c];Counts",60,-1.0,0.0);
	h1p_list.push_back(h2p_prelzPr);
	TH1D * h2p_prelTPr = new TH1D("epp_prelTPr","ep;p_relTPr [GeV/c];Counts",60,0.0,0.5);
	h1p_list.push_back(h2p_prelTPr);
	TH1D * h2p_prel = new TH1D("epp_prel","ep;p_rel [GeV/c];Counts",60,0.0,1.0);
	h1p_list.push_back(h2p_prel);

	TH1D * h2p_e1 =  new TH1D("epp_e1" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1);
	TH1D * h2p_e1_lo =  new TH1D("epp_e1_lo" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1_lo);
	TH1D * h2p_e1_md =  new TH1D("epp_e1_md" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1_md);
	TH1D * h2p_e1_hi =  new TH1D("epp_e1_hi" ,"epp;e1 [GeV];Counts",28,0.3,1.1);
	h1p_list.push_back(h2p_e1_hi);

	
	TH2D * thetaq_vs_phiq = new TH2D("thetaq_vs_phiq","ep; theta_q [degrees]; phi_q [degrees]",40,30.,70.,60,-30.,330.);
	h1p_list.push_back(thetaq_vs_phiq);
	TH2D * theta1_vs_phi1 = new TH2D("theta1_vs_phi1","ep; theta_1 [degrees]; phi_1 [degrees]",80,10.,90.,60,-30.,330.);
	h1p_list.push_back(theta1_vs_phi1);
	TH2D * thetae_vs_phie = new TH2D("thetae_vs_phie","ep; theta_e [degrees]; phi_e [degrees]",60,10.,40.,60,-30.,330.);
	h1p_list.push_back(thetae_vs_phie);

	TH2D * prel_vs_Qsq = new TH2D("prel_vs_Qsq","ep; p_rel [GeV/c]; QSq [GeV^2]",40,0.0,1.0,40,1.,5.);
	h1p_list.push_back(prel_vs_Qsq);
	TH2D * prel_vs_xB = new TH2D("prel_vs_xB","ep; p_rel [GeV/c]; xB",40,0.0,1.0,26,1.2,2.5);
	h1p_list.push_back(prel_vs_xB);
	TH2D * prel_vs_Pm = new TH2D("prel_vs_Pm","ep; p_rel [GeV/c]; Pm [GeV/c]",40,0.0,1.0,28,0.3,1.0);
	h1p_list.push_back(prel_vs_Pm);
	TH2D * prel_vs_Em = new TH2D("prel_vs_Em","ep; p_rel [GeV/c]; Em [GeV]",40,0.0,1.0,40,-0.2,0.6);
	h1p_list.push_back(prel_vs_Em);
	TH2D * prel_vs_cPmq = new TH2D("prel_vs_cPmq","ep; p_rel [GeV/c]; cos(Theta_Pmq)",40,0.0,1.0,40,-1.,0);
	h1p_list.push_back(prel_vs_cPmq);
	TH2D * prel_vs_cPm2 = new TH2D("prel_vs_cPm2","ep; p_rel [GeV/c]; cos(Theta_Pm2)",40,0.0,1.0,40,-1.,0.);
	h1p_list.push_back(prel_vs_cPm2);
	TH2D * prel_vs_mom2 = new TH2D("prel_vs_mom2","ep; p_rel [GeV/c]; mom2 [GeV/c]",40,0.0,1.0,17,0.35,1.2);
	h1p_list.push_back(prel_vs_mom2);

	TH2D * prec_vs_e1 = new TH2D("prec_vs_e1","ep; p_rec [GeV/c]; e1 [GeV]",17,0.35,1.2,40,0.6,1.6);
	h1p_list.push_back(prec_vs_e1);

	TH2D * Pm_vs_Thetam = new TH2D("Pm_vs_Thetam","ep; p_miss [GeV/c]; Theta_m [degrees]",28,0.3,1.0,40,10.,170.);
	h1p_list.push_back(Pm_vs_Thetam);
	TH2D * Pmz_vs_PmT = new TH2D("Pmz_vs_PmT","ep; p_mz [GeV/c]; p_mT [GeV/c]",30,0.,0.4,30,0.0,1.0);
	h1p_list.push_back(Pmz_vs_PmT);
	TH2D * Pmzq_vs_PmTq = new TH2D("Pmzq_vs_PmTq","ep; p_mzq [GeV/c]; p_mTq [GeV/c]",30,-1.0,0.0,30,0.0,1.0);
	h1p_list.push_back(Pmzq_vs_PmTq);

	TH2D * P1_vs_Theta1 = new TH2D("Pm_vs_Theta1","ep; p_1 [GeV/c]; Theta_1 [degrees]",40,0.4,2.4,40,10.,90.);
	h1p_list.push_back(P1_vs_Theta1);
	TH2D * P1z_vs_P1T = new TH2D("P1z_vs_P1T","ep; p_1z [GeV/c]; p_1T [GeV/c]",60,0.0,2.0,60,0.0,2.0);
	h1p_list.push_back(P1z_vs_P1T);
	TH2D * P1zq_vs_P1Tq = new TH2D("P1zq_vs_P1Tq","ep; p_1zq [GeV/c]; p_1Tq [GeV/c]",60,0.5,3.0,60,0.0,1.0);
	h1p_list.push_back(P1zq_vs_P1Tq);

	TH1D * h1p_thetae_bySec[6];
	TH1D * h1p_theta1_bySec[6];
	TH1D * h2p_theta1_bySec[6];
	TH1D * h1p_thetaq_bySec[6];
	TH1D * h1p_mome_bySec[6];
	TH1D * h1p_q_bySece[6];
	TH1D * h1p_mom1_bySec[6];
	TH1D * h1p_cPmq_bySece[6];
	TH1D * h1p_cPmq_bySec[6];
	TH1D * h1p_Pm_bySec[6];
	TH1D * h1p_Pm_bySece[6];
	TH1D * h1p_Em_bySece[6];
	TH1D * h1p_Pm_bySece_tot_10[6];
	TH1D * h1p_Pm_bySece_tot_15[6];
	TH1D * h1p_Pm_bySece_tot_17[6];
	TH1D * h1p_Pm_bySece_tot_20[6];

	for (int i=0 ; i<6 ; i++)
	{
		char temp[100];

		sprintf(temp,"ep_theta1_%d",i);
		h1p_theta1_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h1p_theta1_bySec[i]);

		sprintf(temp,"epp_theta1_%d",i);
		h2p_theta1_bySec[i] = new TH1D(temp,"epp;Theta [deg];Counts",60,10.,130.);
		h1p_list.push_back(h2p_theta1_bySec[i]);

		sprintf(temp,"ep_thetae_%d",i);
		h1p_thetae_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",60,10.,40.);
		h1p_list.push_back(h1p_thetae_bySec[i]);

		sprintf(temp,"ep_thetaq_%d",i);
		h1p_thetaq_bySec[i] = new TH1D(temp,"ep;Theta [deg];Counts",40,30.,70.);
		h1p_list.push_back(h1p_thetaq_bySec[i]);

		sprintf(temp,"ep_mome_%d",i);
		h1p_mome_bySec[i] = new TH1D(temp,"ep;Mom_e [GeV/c];Counts",40,3.0,5.0);
		h1p_list.push_back(h1p_mome_bySec[i]);

		sprintf(temp,"ep_mom1_%d",i);
		h1p_mom1_bySec[i] = new TH1D(temp,"ep;Mom_1 [GeV/c];Counts",40,0.4,2.4);
		h1p_list.push_back(h1p_mom1_bySec[i]);

		sprintf(temp,"ep_q_%d",i);
		h1p_q_bySece[i] = new TH1D(temp,"ep;q [GeV/c];Counts",40,1.0,3.0);
		h1p_list.push_back(h1p_q_bySece[i]);

		sprintf(temp,"ep_cPmq_e_%d",i);
		h1p_cPmq_bySece[i] = new TH1D(temp,"ep;cos(Theta_Pmq);Counts",40,-1.,0.);
		h1p_list.push_back(h1p_cPmq_bySece[i]);

		sprintf(temp,"ep_cPmq_%d",i);
		h1p_cPmq_bySec[i] = new TH1D(temp,"ep;cos(Theta_Pmq);Counts",40,-1.,0.);
		h1p_list.push_back(h1p_cPmq_bySec[i]);

		sprintf(temp,"ep_Pm_%d",i);
		h1p_Pm_bySec[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySec[i]);

		sprintf(temp,"ep_Pm_e_%d",i);
		h1p_Pm_bySece[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece[i]);

		sprintf(temp,"ep_Em_e_%d",i);
		h1p_Em_bySece[i] = new TH1D(temp,"ep;EMiss [GeV];Counts",40,-0.2,0.6);
		h1p_list.push_back(h1p_Em_bySece[i]);

		sprintf(temp,"ep_Pm_etot_10_%d",i);
		h1p_Pm_bySece_tot_10[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_10[i]);

		sprintf(temp,"ep_Pm_etot_15_%d",i);
		h1p_Pm_bySece_tot_15[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_15[i]);

		sprintf(temp,"ep_Pm_etot_17_%d",i);
		h1p_Pm_bySece_tot_17[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_17[i]);

		sprintf(temp,"ep_Pm_etot_20_%d",i);
		h1p_Pm_bySece_tot_20[i] = new TH1D(temp,"ep;pMiss [GeV];Counts",28,0.3,1.0);
		h1p_list.push_back(h1p_Pm_bySece_tot_20[i]);

	}

	// Now that all histograms have been defined, set them to Sumw2
	for (int i=0 ; i<h1p_list.size() ; i++)
	  h1p_list[i]->Sumw2();

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
		if (Pmiss_size[0] < pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;

		// Let's make a sanitized phi and sector
		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vp.Theta() * 180./M_PI;

		if ((sector==90) and (theta1_deg < 52))
		  continue;
		if ((sector==92) and ((theta1_deg < 45) or (theta1_deg > 72)))
		  continue;
		if ((sector==93) and (theta1_deg < 50))
		  continue;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);

		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = clas_sector(phie_deg);
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 qu = vq.Unit();

		double phiq_deg = vq.Phi() * 180./M_PI;
		if (phiq_deg < -30.)
		  phiq_deg += 360.;
		double thetaq_deg = vq.Theta() * 180./M_PI;

		h1p_thetaq->Fill(thetaq_deg,weight);
		h1p_q->Fill(vq.Mag(),weight);
		h1p_qz ->Fill(vq.z(),weight);
		h1p_qT ->Fill(vq.XYvector().Mod(),weight);

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		double omega = Q2/(2.*mN*Xb);

		TVector3 vm = vp - vq;
		double thetam_deg = vm.Theta() * 180./M_PI;

		h1p_thetam->Fill(thetam_deg,weight);
		h1p_cthetam->Fill(cos(thetam_deg * M_PI/180.),weight);

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);

		h1p_mome->Fill(ve.Mag(),weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		h1p_momez ->Fill(ve.z(),weight);
		h1p_momeT ->Fill(ve.XYvector().Mod(),weight);
		h1p_mom1z ->Fill(vp.z(),weight);
		h1p_mom1T ->Fill(vp.XYvector().Mod(),weight);
		h1p_Pmz ->Fill(vm.z(),weight);
		h1p_PmT ->Fill(vm.XYvector().Mod(),weight);

		h1p_momezq ->Fill(ve.Dot(qu),weight);
		h1p_momeTq ->Fill(ve.Perp(qu),weight);
		h1p_mom1zq ->Fill(vp.Dot(qu),weight);
		h1p_mom1Tq ->Fill(vp.Perp(qu),weight);
		h1p_Pmzq ->Fill(vm.Dot(qu),weight);
		h1p_PmTq ->Fill(vm.Perp(qu),weight);

		h1p_cPeq->Fill(cos(ve.Angle(vq)),weight);
		h1p_cP1q->Fill(cos(vp.Angle(vq)),weight);
		h1p_cPm1->Fill(cos(vm.Angle(vp)),weight);

		Pm_vs_Thetam->Fill(Pmiss_size[0],thetam_deg,weight);
		Pmz_vs_PmT->Fill(vm.z(),vm.XYvector().Mod(),weight);
		Pmzq_vs_PmTq->Fill(vm.Dot(qu),vm.Perp(qu),weight);
		
		P1_vs_Theta1->Fill(Pp_size[0],theta1_deg,weight);
		P1z_vs_P1T->Fill(vp.z(),vp.XYvector().Mod(),weight);
		P1zq_vs_P1Tq->Fill(vp.Dot(qu),vp.Perp(qu),weight);

		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		h1p_thetaq_bySec[sec_e]->Fill(thetaq_deg,weight);
		h1p_q_bySece[sec_e]->Fill(vq.Mag(),weight);
		h1p_cPmq_bySece[sec_e]->Fill(cos(Pmiss_q_angle[0]*M_PI/180),weight); 
		h1p_cPmq_bySec[sector]->Fill(cos(Pmiss_q_angle[0]*M_PI/180),weight); 
		h1p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		h1p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 
		h1p_Em_bySece[sec_e]->Fill(Emiss,weight); 		
		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1_bySec[sector]->Fill(Pp_size[0],weight);

		thetaq_vs_phiq->Fill(thetaq_deg,phiq_deg,weight);
		theta1_vs_phi1->Fill(theta1_deg,phi1_deg,weight);
		thetae_vs_phie->Fill(thetae_deg,phie_deg,weight);


		if (thetae_deg > 10)
		  {
		    h1p_Pm_10 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_10[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 15)
		  {
		    h1p_Pm_15 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_15[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 17.5)
		  {
		    h1p_Pm_17 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_17[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 20)
		  {
		    h1p_Pm_20 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_20[sec_e] ->Fill(Pmiss_size[0],weight);
		  }

	}

	f1p->Close();

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
		if (Pmiss_size[0] < pmiss_cut)
		  continue;

		// Apply fiducial cuts
		TVector3 ve(Pe[0],Pe[1],Pe[2]);
		TVector3 vp(Pp[0][0],Pp[0][1],Pp[0][2]);
		if (!accept_electron(ve))
		  continue;
		if (!accept_proton(vp))
		  continue;

		// Let's make a sanitized phi and sector
		double phi1_deg = vp.Phi() * 180./M_PI;
		if (phi1_deg < -30.)
			phi1_deg += 360.;
		int sector = clas_sector(phi1_deg);
		double theta1_deg = vp.Theta() * 180./M_PI;

		if ((sector==90) and (theta1_deg < 52))
		  continue;
		if ((sector==92) and ((theta1_deg < 45) or (theta1_deg > 72)))
		  continue;
		if ((sector==93) and (theta1_deg < 50))
		  continue;

		h1p_phi1->Fill(phi1_deg,weight);
		h1p_theta1->Fill(theta1_deg,weight);

		double phie_deg = ve.Phi() * 180./M_PI;
		if (phie_deg < -30.)
			phie_deg += 360.;
		int sec_e = clas_sector(phie_deg);
		double thetae_deg = ve.Theta() * 180./M_PI;

		h1p_phie->Fill(phie_deg,weight);
		h1p_thetae->Fill(thetae_deg,weight);

		TVector3 vq(q[0],q[1],q[2]);
		TVector3 qu = vq.Unit();

		double phiq_deg = vq.Phi() * 180./M_PI;
		if (phiq_deg < -30.)
		  phiq_deg += 360.;
		double thetaq_deg = vq.Theta() * 180./M_PI;

		h1p_thetaq->Fill(thetaq_deg,weight);
		h1p_q->Fill(vq.Mag(),weight);
		h1p_qz ->Fill(vq.z(),weight);
		h1p_qT ->Fill(vq.XYvector().Mod(),weight);

		h1p_QSq->Fill(Q2,weight);
		h1p_xB ->Fill(Xb,weight);
		h1p_Pm ->Fill(Pmiss_size[0],weight);
		h1p_cPmq->Fill(cos(Pmiss_q_angle[0]*M_PI/180.),weight);

		double omega = Q2/(2.*mN*Xb);
		double e1 = sqrt(vp.Mag2()+mN*mN)-omega;

		TVector3 vm = vp - vq;
		TVector3 mu = vm.Unit();
		double thetam_deg = vm.Theta() * 180./M_PI;

		h1p_thetam->Fill(thetam_deg,weight);
		h1p_cthetam->Fill(cos(thetam_deg * M_PI/180.),weight);

		double Emiss = Q2/(2.*mN*Xb) + m_12C - sqrt(Pp_size[0]*Pp_size[0] + mN*mN) - sqrt(Pmiss_size[0]*Pmiss_size[0] + m_11B*m_11B);

		h1p_mome->Fill(ve.Mag(),weight);
		h1p_mom1->Fill(Pp_size[0],weight);
		h1p_momez ->Fill(ve.z(),weight);
		h1p_momeT ->Fill(ve.XYvector().Mod(),weight);
		h1p_mom1z ->Fill(vp.z(),weight);
		h1p_mom1T ->Fill(vp.XYvector().Mod(),weight);
		h1p_Pmz ->Fill(vm.z(),weight);
		h1p_PmT ->Fill(vm.XYvector().Mod(),weight);

		h1p_momezq ->Fill(ve.Dot(qu),weight);
		h1p_momeTq ->Fill(ve.Perp(qu),weight);
		h1p_mom1zq ->Fill(vp.Dot(qu),weight);
		h1p_mom1Tq ->Fill(vp.Perp(qu),weight);
		h1p_Pmzq ->Fill(vm.Dot(qu),weight);
		h1p_PmTq ->Fill(vm.Perp(qu),weight);

		h1p_cPeq->Fill(cos(ve.Angle(vq)),weight);
		h1p_cP1q->Fill(cos(vp.Angle(vq)),weight);
		h1p_cPm1->Fill(cos(vm.Angle(vp)),weight);

		Pm_vs_Thetam->Fill(Pmiss_size[0],thetam_deg,weight);
		Pmz_vs_PmT->Fill(vm.z(),vm.XYvector().Mod(),weight);
		Pmzq_vs_PmTq->Fill(vm.Dot(qu),vm.Perp(qu),weight);
		
		P1_vs_Theta1->Fill(Pp_size[0],theta1_deg,weight);
		P1z_vs_P1T->Fill(vp.z(),vp.XYvector().Mod(),weight);
		P1zq_vs_P1Tq->Fill(vp.Dot(qu),vp.Perp(qu),weight);

		h1p_thetae_bySec[sec_e]->Fill(thetae_deg,weight);
		h1p_mome_bySec[sec_e]->Fill(ve.Mag(),weight);
		h1p_thetaq_bySec[sec_e]->Fill(thetaq_deg,weight);
		h1p_q_bySece[sec_e]->Fill(vq.Mag(),weight);
		h1p_cPmq_bySece[sec_e]->Fill(cos(Pmiss_q_angle[0]*M_PI/180),weight); 
		h1p_cPmq_bySec[sec_e]->Fill(cos(Pmiss_q_angle[0]*M_PI/180),weight); 
		h1p_Pm_bySec[sector]->Fill(Pmiss_size[0],weight); 
		h1p_Pm_bySece[sec_e]->Fill(Pmiss_size[0],weight); 
		h1p_Em_bySece[sec_e]->Fill(Emiss,weight); 

		h1p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h1p_mom1_bySec[sector]->Fill(Pp_size[0],weight);

		thetaq_vs_phiq->Fill(thetaq_deg,phiq_deg,weight);
		theta1_vs_phi1->Fill(theta1_deg,phi1_deg,weight);
		thetae_vs_phie->Fill(thetae_deg,phie_deg,weight);

		if (thetae_deg > 10)
		  {
		    h1p_Pm_10 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_10[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 15)
		  {
		    h1p_Pm_15 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_15[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 17.5)
		  {
		    h1p_Pm_17 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_17[sec_e] ->Fill(Pmiss_size[0],weight);
		  }
		if (thetae_deg > 20)
		  {
		    h1p_Pm_20 ->Fill(Pmiss_size[0],weight);
		    h1p_Pm_bySece_tot_20[sec_e] ->Fill(Pmiss_size[0],weight);
		  }

		if (Pp_size[1] < 0.35)
		  continue;
		h2p_Pm ->Fill(Pmiss_size[0],weight);

		h2p_e1 ->Fill(e1,weight);
		if (Pp_size[1] < 0.45)
		  h2p_e1_lo ->Fill(e1,weight);
		else if (Pp_size[1] < 0.55)
		  h2p_e1_md ->Fill(e1,weight);
		else
		  h2p_e1_hi ->Fill(e1,weight);

		h2p_theta1_bySec[sector]->Fill(theta1_deg,weight);
		h2p_theta1->Fill(theta1_deg,weight);

		TVector3 vrec(Pp[1][0],Pp[1][1],Pp[1][2]);
		TVector3 ru = vrec.Unit();
		h2p_mom2z ->Fill(vrec.z(),weight);
		h2p_mom2T ->Fill(vrec.XYvector().Mod(),weight);
		h2p_mom2zq ->Fill(vrec.Dot(qu),weight);
		h2p_mom2Tq ->Fill(vrec.Perp(qu),weight);
		
		TVector3 vcm = vm + vrec;
		TVector3 ucm = vcm.Unit();
		TVector3 vrel = (vm - vrec)*(0.5);
		TVector3 urel = vrel.Unit();
		Double_t prel = vrel.Mag();
		h2p_pCM ->Fill(vcm.Mag(),weight);
		h2p_pCMz ->Fill(vcm.z(),weight);
		h2p_pCMT ->Fill(vcm.XYvector().Mod(),weight);
		h2p_pCMzPm ->Fill(vcm.Dot(mu),weight);
		h2p_pCMTPm ->Fill(vcm.Perp(mu),weight);
		h2p_pCMzq ->Fill(vcm.Dot(qu),weight);
		h2p_pCMTq ->Fill(vcm.Perp(qu),weight);
		h2p_pCMzPr ->Fill(vcm.Dot(ru),weight);
		h2p_pCMTPr ->Fill(vcm.Perp(ru),weight);
		h2p_prel ->Fill(prel,weight);
		h2p_prelz ->Fill(vrel.z(),weight);
		h2p_prelT ->Fill(vrel.XYvector().Mod(),weight);
		h2p_prelzPm ->Fill(vrel.Dot(mu),weight);
		h2p_prelTPm ->Fill(vrel.Perp(mu),weight);
		h2p_prelzq ->Fill(vrel.Dot(qu),weight);
		h2p_prelTq ->Fill(vrel.Perp(qu),weight);
		h2p_prelzPr ->Fill(vrel.Dot(ru),weight);
		h2p_prelTPr ->Fill(vrel.Perp(ru),weight);

		h2p_Pmzrel ->Fill(vm.Dot(urel),weight);
		h2p_PmTrel ->Fill(vm.Perp(urel),weight);
		h2p_Pmzrec ->Fill(vm.Dot(ru),weight);
		h2p_PmTrec ->Fill(vm.Perp(ru),weight);
		h2p_PmzCM ->Fill(vm.Dot(ucm),weight);
		h2p_PmTCM ->Fill(vm.Perp(ucm),weight);
		
		h2p_cPm2 ->Fill(cos(vm.Angle(vrec)),weight);

		prel_vs_Qsq->Fill(prel,Q2,weight);
		prel_vs_xB->Fill(prel,Xb,weight);
		prel_vs_Pm->Fill(prel,Pmiss_size[0],weight);
		prel_vs_Em->Fill(prel,Emiss,weight);
		prel_vs_cPmq->Fill(prel,cos(Pmiss_q_angle[0]*M_PI/180),weight);
		prel_vs_cPm2->Fill(prel,cos(vm.Angle(vrec)),weight);
		prel_vs_mom2->Fill(prel,Pp_size[1],weight);

		prec_vs_e1->Fill(Pp_size[1],e1,weight);

		h2p_mom2->Fill(Pp_size[1],weight);

	}

	f2p->Close();

	cerr << "The ep integrals is: " << h1p_Pm->Integral() << "\n";
	cerr << "Broken down by pmiss range...\n\n";
	for (int j=4 ; j<=28 ; j+=4)
	  {
	    double min=0.3 + 0.1*(j-4)/4.;
	    double max=0.3 + 0.1*j/4.;
	    cerr << min << " < pmiss < " << max << " : " << h1p_Pm->Integral(j-3,j) << "\n";
	  }

	// Write out
	fo -> cd();

	const double data_ep = 4945.;
	const double data_epp = 359.;
	const double pnorm = data_ep/h1p_Pm->Integral();

	// scale all the histograms, and write them out
	for (int i=0 ; i<h1p_list.size() ; i++)
	  {
	    h1p_list[i]->Scale(pnorm);
	    h1p_list[i]->Write();
	  }

	fo->Close();

	return 0;
}

#ifndef __NUCLEAR_INFO_H__
#define __NUCLEAR_INFO_H__
#include <cstdlib>
#include <vector>

// EG2 info
const double eg2beam=5.014;

// constants
const double me = 0.000511;
const double mN = 0.93892;
const double mU=0.931;
const double GeVfm=0.1973;
const double alpha=0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;

// nuclear masses
const double m_1H = mN;
const double m_2H = 2.01410178*mU;
const double m_4He = 4.00260325415 * mU;
const double m_6Li = 6.015122795 * mU;
const double m_8Be = 8.00530510 * mU;
const double m_10B = 10.0129370 * mU;
const double m_11B = 11.0093054 * mU;
const double m_12C = 12.*mU;
const double m_14N = 14.0030740048*mU;
const double m_16O = 15.99491461956*mU;

// nucleon codes
const int pCode=2212;
const int nCode=2112;

class Nuclear_Info
{
 public:
  Nuclear_Info(int thisA, int pType);
  ~Nuclear_Info();
  void set_Estar(double new_Estar);
  void set_Contacts(double new_Cpp0, double new_Cpn0, double new_Cpn1);
  
  double get_Estar();
  double get_Cpp0();
  double get_Cpn0();
  double get_Cpn1();
  
  double get_pp(double k_rel);
  double get_pn(double k_rel);
  double get_pn0(double k_rel);
  double get_pn1(double k_rel);
  double get_mA();
  double get_mAm2();
  double get_sigmaCM();
  void do_SXC(int &lead_type, int &rec_type, double r);
  // std::vector<double> get_SCX_Ps();

 private:
  int A;
  double mA;
  double mAm2;
  double Estar;// = 0;
  double sigmaCM;
  double d_sigmaCM;// = 0.;
  double phiSq_pp0[100];
  double phiSq_pn0[100];
  double phiSq_pn1[100];
  
  double Cpp0;
  double d_Cpp0;// = 0.;
  double Cpn0;
  double d_Cpn0;// = 0.;
  double Cpn1;
  double d_Cpn1;// = 0.;
  
  double get_phiSq(double *phiPtr, double k_rel);
  void fill_arrays();
  void fill_arrays_chiral();
  void fill_arrays_chiral_n3lo();

  double pPP2PN;// = 0.;
  double d_pPP2PN;// = 0.;
  double pPP2NP;// = 0.;
  double d_pPP2NP;// = 0.;
  double pPP2NN;// = 0.;
  double d_pPP2NN;// = 0.;

  double pPN2NN;// = 0.;
  double d_pPN2NN;// = 0;
  double pPN2PP;// = 0.;
  double d_pPN2PP;// = 0.;
  double pPN2NP;// = 0.;
  double d_pPN2NP;// = 0.;

  double pNP2NN;// = 0.;
  double d_pNP2NN;// = 0.;
  double pNP2PP;// = 0.;
  double d_pNP2PP;// = 0.;
  double pNP2PN;// = 0.;
  double d_pNP2PN;// = 0.;

  double pNN2PN;// = 0.;
  double d_pNN2PN;// = 0.;
  double pNN2NP;// = 0.;
  double d_pNN2NP;// = 0.;
  double pNN2PP;// = 0.;
  double d_pNN2PP;// = 0.;


};

#endif

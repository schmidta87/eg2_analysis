 #ifndef __NUCLEAR_INFO_H__
 #define __NUCLEAR_INFO_H__
 #include <cstdlib>
 #include <vector>
 #include "constants.h"

 // nuclear masses
 const double m_1H = mN;
 const double m_2H = 2.01410178*mU - me;
 const double m_4He = 4.00260325415 * mU - 2*me;
 const double m_6Li = 6.015122795 * mU - 3*me;
 const double m_8Be = 8.00530510 * mU - 4*me;
 const double m_10B = 10.0129370 * mU - 5*me;
 const double m_11B = 11.0093054 * mU - 5*me;
 const double m_12C = 12.*mU - 6*me;
 const double m_14N = 14.0030740048*mU - 7*me;
 const double m_16O = 15.99491461956*mU - 8*me;

 // nucleon codes
 const int pCode=2212;
 const int nCode=2112;

 class Nuclear_Info
 {
  public:
   Nuclear_Info(int thisA, int pType);
   ~Nuclear_Info();
   void set_sigmaCM(double new_sigmaCM);
   void set_Estar(double new_Estar);
   void set_Contacts(double new_Cpp0, double new_Cpn0, double new_Cpn1);
   void set_byRatio();
   void set_Ratio(double ratio);
   void set_Ratio(double ratio, double sig_ratio);
   void set_Ratio_Inverse(double ratio);
   void set_Ratio_Inverse(double ratio, double sig_ratio);
   void randomize();
   
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
   std::vector<double> get_SCX_Ps();

 private:
  int A;
  int u;
  double mA;
  double mAm2;
  double Estar;
  double sigmaCM;
  double d_sigmaCM;
  double phiSq_pp0[100];
  double phiSq_pn0[100];
  double phiSq_pn1[100];
  
  double Cpp0;
  double d_Cpp0;
  double Cpn0;
  double d_Cpn0;
  double Cpn1;
  double d_Cpn1;
  
  double get_phiSq(double *phiPtr, double k_rel);
  void fill_arrays();
  void fill_arrays_chiral();
  void fill_arrays_chiral_n3lo();
  void fill_arrays_chiral_R12();

  double pPP2NP;
  double d_pPP2NP;
  double pPP2PN;
  double d_pPP2PN;
  double pPP2NN;
  double d_pPP2NN;

  double pPN2NN;
  double d_pPN2NN;
  double pPN2PP;
  double d_pPN2PP;
  double pPN2NP;
  double d_pPN2NP;

  double pNP2PP;
  double d_pNP2PP;
  double pNP2NN;
  double d_pNP2NN;
  double pNP2PN;
  double d_pNP2PN;

  double pNN2PN;
  double d_pNN2PN;
  double pNN2NP;
  double d_pNN2NP;
  double pNN2PP;
  double d_pNN2PP;

};

#endif

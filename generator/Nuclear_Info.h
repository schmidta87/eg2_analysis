#ifndef __NUCLEAR_INFO_H__
#define __NUCLEAR_INFO_H__

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

class Nuclear_Info
{
 public:
  Nuclear_Info(int thisA);
  ~Nuclear_Info();
  double get_pp(double k_rel);
  double get_pn(double k_rel);
  double get_pn0(double k_rel);
  double get_pn1(double k_rel);
  double get_mA();
  double get_mAm2();
  double get_sigmaCM();

 private:
  int A;
  double mA;
  double mAm2;
  double sigmaCM;
  double phiSq_pp0[100];
  double phiSq_pn0[100];
  double phiSq_pn1[100];
  double Cpp0;
  double Cpn0;
  double Cpn1;
  double get_phiSq(double *phiPtr, double k_rel);
  void fill_arrays();
  void fill_arrays_chiral();
  
};

#endif

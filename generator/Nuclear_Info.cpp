#include "Nuclear_Info.h"

#include "constants.h"

#include <iostream>
#include <cstdlib>
#include <vector>

Nuclear_Info::Nuclear_Info(int thisA, int pType)
{
  A = thisA;

  if(pType == 1){
  fill_arrays();
  std::cerr <<"You are using the AV18 potential\n";
  }
  else if (pType == 2){
    fill_arrays_chiral();
    std::cerr <<"You are using the N2L0 potential\n";
  }
  else if (pType == 3){
    fill_arrays_chiral_n3lo();
    std::cerr <<"You are using the N3L0 potential\n";
    
  }
  else{
    std::cerr <<"You are using a potential not in the library. \n Aborting...\n";
  exit(-2);
  }

  Estar = 0;

  if (A==2)
    {
      sigmaCM=0.;
      Cpp0=0.;
      Cpn0=0.;
      Cpn1=4.;
      mA=m_2H;
      mAm2=0.;
    }
  else if (A==4)
    {
      //sigmaCM=0.;
      //Cpp0=0.;
      //Cpn0=0.;
      //Cpn1=16.0;
      //mA=m_2H;
      //mAm2=m_2H;

      sigmaCM=0.1;
      Cpp0=0.65;
      Cpn0=0.69;
      Cpn1=12.3;
      mA=m_4He;
      mAm2=m_2H;
    }
  else if (A==12)
    {
      //sigmaCM=0.;
      //Cpp0 = 0.;
      //Cpn0 = 0.;
      //Cpn1 = 16.;
      //mA = m_12C;
      //mAm2 = m_10B;

      sigmaCM=0.15;
      d_sigmaCM = 0.02;
      Cpp0 = 1.3;
      d_Cpp0 = 0.2;
      Cpn0 = 1.4;
      d_Cpn0 = 0.2;
      Cpn1 = 16.8;
      d_Cpn1 = 0.8;
      mA = m_12C;
      mAm2 = m_10B;
      
      pPP2PN = 0.048;
      d_pPP2PN = 0.003;
      pPP2NP = 0.041;
      d_pPP2NP = 0.003;
      pPP2NN = 0.0029;
      d_pPP2NN = 0.0002;
      
      pPN2NN = 0.035;
      d_pPN2NN = 0.002;
      pPN2PP = 0.041;
      d_pPN2PP = 0.003;
      pPN2NP = 0.0021;
      d_pPN2NP = 0.0001;
      
      pNP2NN = 0.041;
      d_pNP2NN = 0.003;
      pNP2PP = 0.035;
      d_pNP2PP = 0.002;
      pNP2PN = 0.0021;
      d_pNP2PN = 0.0001;
      
      pNN2PN = 0.041;
      d_pNN2PN = 0.003;
      pNN2NP = 0.048;
      d_pNN2NP = 0.003;
      pNN2PP = 0.0029;
      d_pNN2PP = 0.0002;
    }
  else
    {
      std::cerr << "You selected a nucleus with A=" << A << "\n"
	   << " which is not in the library. Aborting...\n";
      exit(-2);
    }
}

Nuclear_Info::~Nuclear_Info()
{
}

void Nuclear_Info::set_sigmaCM(double new_sigmaCM)
{
  sigmaCM = new_sigmaCM;
}

void Nuclear_Info::set_Estar(double new_Estar)
{
  Estar = new_Estar;
}

void Nuclear_Info::set_Contacts(double new_Cpp0, double new_Cpn0, double new_Cpn1)
{
  Cpp0 = new_Cpp0;
  Cpn0 = new_Cpn0;
  Cpn1 = new_Cpn1;
}

void Nuclear_Info::randomize()
{
}

double Nuclear_Info::get_Estar()
{
  return Estar;
}

double Nuclear_Info::get_Cpp0()
{
  return Cpp0;
}

double Nuclear_Info::get_Cpn0()
{
  return Cpn0;
}

double Nuclear_Info::get_Cpn1()
{
  return Cpn1;
}

double Nuclear_Info::get_mA()
{
  return mA;
}

double Nuclear_Info::get_mAm2()
{
  return mAm2 + Estar;
}

double Nuclear_Info::get_sigmaCM()
{
  return sigmaCM;
}

double Nuclear_Info::get_pp(double k_rel)
{
  return 2. * 0.5*A*Cpp0 * get_phiSq(phiSq_pp0,k_rel); // The 2 comes from contact definition
}

double Nuclear_Info::get_pn(double k_rel)
{
  return get_pn0(k_rel) + get_pn1(k_rel);
}

double Nuclear_Info::get_pn0(double k_rel)
{
  return 0.5*A*Cpn0 * get_phiSq(phiSq_pn0,k_rel);
}

double Nuclear_Info::get_pn1(double k_rel)
{
  return 0.5*A*Cpn1 * get_phiSq(phiSq_pn1,k_rel);
}

void Nuclear_Info::do_SXC(int &lead_type, int &rec_type, double r)
{
  // Now we do a whole bunch of tests
  if ((lead_type==pCode) && (rec_type==pCode))
    {
      if (r < pPP2NN)
	{
	  lead_type=nCode;
	  rec_type=nCode;
	}
      else if (r < pPP2NN + pPP2NP)
	{
	  lead_type=nCode;
	}
      else if (r < pPP2NN + pPP2NP + pPP2PN)
	{
	  rec_type=nCode;
	}
    }
  else if ((lead_type==nCode) && (rec_type==pCode))
    {
      if (r < pNP2PN)
	{
	  lead_type=pCode;
	  rec_type=nCode;
	}
      else if (r < pNP2PN + pNP2PP)
	{
	  lead_type=pCode;
	}
      else if (r < pNP2PN + pNP2PP + pNP2NN)
	{
	  rec_type=nCode;
	}
    }
  else if ((lead_type==pCode) && (rec_type==nCode))
    {
      if (r < pPN2NP)
	{
	  lead_type=nCode;
	  rec_type=pCode;
	}
      else if (r < pPN2NP + pPN2NN)
	{
	  lead_type=nCode;
	}
      else if (r < pPN2NP + pPN2NN + pPN2PP)
	{
	  rec_type=pCode;
	}
    }
  else if ((lead_type==nCode) && (rec_type==nCode))
    {
      if (r < pNN2PP)
	{
	  lead_type=pCode;
	  rec_type=pCode;
	}
      else if (r < pNN2PP + pNN2PN)
	{
	  lead_type=pCode;
	}
      else if (r < pNN2PP + pNN2PN + pNN2NP)
	{
	  rec_type=pCode;
	}
    }
  else
    {
      std::cerr << "Invalid nucleon codes. Check and fix. Exiting\n\n\n";
    }
}

std::vector<double> Nuclear_Info::get_SCX_Ps()
{
  std::vector<double> Ps;
  Ps.push_back(pPP2PN);
  Ps.push_back(pPP2NP);
  Ps.push_back(pPP2NN);
  Ps.push_back(pPN2NN);
  Ps.push_back(pPN2PP);
  Ps.push_back(pPN2NP);
  Ps.push_back(pNP2NN);
  Ps.push_back(pNP2PP);
  Ps.push_back(pNP2PN);
  Ps.push_back(pNN2PN);
  Ps.push_back(pNN2NP);
  Ps.push_back(pNN2PP);
  return Ps;
}

double Nuclear_Info::get_phiSq(double *phiPtr, double k_rel)
{
  double bin = k_rel / GeVfm / 0.1;

  if (bin < 0.)
    return 0.;
  if (bin < 1.)
    return bin * phiPtr[0];
  if (bin > 100.)
    return 0.;
  
  int b = bin;
  double x = bin - b;
  return x*phiPtr[b] + (1.-x)*phiPtr[b-1];
}

void Nuclear_Info::fill_arrays()
{
  double temp_pp0[100]={
18683396.6234051,
639762.771794429,
73908.4475109856,
19752.8915725184,
6985.74841560849,
2883.66695171470,
1311.78054811353,
635.881495153820,
321.191776313180,
166.247451460222,
86.9438198909215,
45.3289533225080,
23.2131420726046,
11.4576914568038,
5.29974612988386,
2.18689765695214,
0.723732923792566,
0.137896096570853,
4.12685844627030e-05,
0.0727158808737714,
0.226901073028930,
0.395063240611470,
0.544381242558901,
0.661344392334904,
0.742861516856787,
0.791143734434978,
0.810792168659347,
0.807178799717402,
0.785584610798292,
0.750777638541285,
0.706842253126365,
0.657147414337408,
0.604387561560645,
0.550657368958016,
0.497538209703937,
0.446184204150935,
0.397401619319003,
0.351718915091108,
0.309446692517894,
0.270727876691614,
0.235579021084050,
0.203923808405533,
0.175619887900554,
0.150480111027402,
0.128289127363092,
0.108816186511540,
0.0918248549792178,
0.0770802581221853,
0.0643543427291766,
0.0534295732034895,
0.0441013995568204,
0.0361797665632472,
0.0294898887833488,
0.0238724658630861,
0.0191834816190184,
0.0152937004965316,
0.0120879498074588,
0.00946426056443420,
0.00733292073493141,
0.00561548499534776,
0.00424377382159058,
0.00315888638862199,
0.00231024647415919,
0.00165469372492237,
0.00115562981773023,
0.000782224904061435,
0.000508687260620702,
0.000313597316489482,
0.000179305093435832,
9.13894307095126e-05,
3.81761072994801e-05,
1.03115213746831e-05,
3.88175399680274e-07,
2.61789357516966e-06,
1.25486971126193e-05,
2.68211999506852e-05,
4.29605262757271e-05,
5.91998708366808e-05,
7.43320801459427e-05,
8.75857955193904e-05,
9.85230339935429e-05,
0.000106955293396467,
0.000112875544007094,
0.000116403787155677,
0.000117744015233489,
0.000117150772586753,
0.000114903644185344,
0.000111288261107098,
0.000106582623299399,
0.000101047641953342,
9.49210697898548e-05,
8.84140257395268e-05,
8.17095142128591e-05,
7.49624331697077e-05,
6.83006230706950e-05,
6.18266636168423e-05,
5.56201120806382e-05,
4.97400019576200e-05,
4.42274383700330e-05,
3.91081590235919e-05 };
  
  double temp_pn0[100]={
    15300347.6134154,
605303.564201820,
73834.3503740097,
19878.2097394933,
7079.24825845492,
2944.09198372589,
1349.97657879891,
660.122202477167,
336.728804141499,
176.298704218147,
93.4845047260791,
49.5886951870085,
25.9714168385921,
13.2178860083239,
6.39260311945783,
2.83299355960728,
1.07195373878876,
0.289625746873957,
0.0246844530122985,
0.0178583062118798,
0.125082288431188,
0.268366005728782,
0.407585843596093,
0.524122463013246,
0.611316245840212,
0.668898684124935,
0.699764403042386,
0.708128892239602,
0.698507640769719,
0.675180174807435,
0.641937252769538,
0.601989917922953,
0.557967692823637,
0.511962616755094,
0.465593717446985,
0.420077426590986,
0.376295990718513,
0.334859908218734,
0.296162721663445,
0.260427816755645,
0.227747609785500,
0.198115831807818,
0.171453776490820,
0.147631376716036,
0.126483928292068,
0.107825202141938,
0.0914575815220566,
0.0771797811605364,
0.0647926089294318,
0.0541031584298509,
0.0449277547493419,
0.0370939135546064,
0.0304415320756220,
0.0248234844428549,
0.0201057640952047,
0.0161672877188262,
0.0128994510382827,
0.0102055111519582,
0.00799985200444371,
0.00620717943986927,
0.00476168126061856,
0.00360617928917502,
0.00269129482205885,
0.00197464213155372,
0.00142006146748792,
0.000996898813714605,
0.000679337007646729,
0.000445780871960240,
0.000278296795707633,
0.000162106335125849,
8.51320544200502e-05,
3.75932352004908e-05,
1.16485714068543e-05,
1.08253163070344e-06,
1.03199417632785e-06,
7.74959578942479e-06,
1.84003059687360e-05,
3.08877976502499e-05,
4.37073532356883e-05,
5.58221914364596e-05,
6.65603459094842e-05,
7.55294202452235e-05,
8.25467753782413e-05,
8.75829716572109e-05,
9.07164467431434e-05,
9.20977102440348e-05,
9.19214729919702e-05,
9.04053578688611e-05,
8.77740288782096e-05,
8.42476860300348e-05,
8.00340987383894e-05,
7.53234078245167e-05,
7.02850949857751e-05,
6.50666109491192e-05,
5.97932169009648e-05,
5.45687270792788e-05,
4.94768459445548e-05,
4.45829003113426e-05,
3.99357918975429e-05,
3.55700295210602e-05 };

 double temp_pn1[100]={
9719347.01980097,
318593.292678900,
22749.3998006928,
6401.62362578146,
2348.52221536720,
1009.39011304670,
480.992840533724,
246.667535493963,
133.815565263115,
76.0204380885667,
44.9698215549149,
27.6256104398493,
17.6113651545925,
11.6554909395860,
8.01265456085003,
5.72141474729729,
4.23766865489032,
3.24648817294686,
2.56196536296042,
2.07249436009417,
1.71001251328634,
1.43240940954246,
1.21323306249160,
1.03560664263678,
0.888560402592821,
0.764807843162692,
0.659386766232806,
0.568808616443671,
0.490541162045559,
0.422675899627828,
0.363724250528384,
0.312484667341479,
0.267962662929490,
0.229313012167132,
0.195808646790267,
0.166814491586681,
0.141772688832050,
0.120190610300630,
0.101632335812215,
0.0857119717105394,
0.0720878567712070,
0.0604582089317134,
0.0505566821510119,
0.0421488714416281,
0.0350288866094745,
0.0290163933649339,
0.0239537816252370,
0.0197037456548734,
0.0161468712195542,
0.0131797264665244,
0.0107128045439877,
0.00866898998945952,
0.00698195608662072,
0.00559481879189304,
0.00445896436216796,
0.00353291789332953,
0.00278145326361891,
0.00217468992210455,
0.00168740683628140,
0.00129835696471829,
0.000989711865063592,
0.000746569317164056,
0.000556510241080911,
0.000409234833246501,
0.000296231028763088,
0.000210498717808956,
0.000146307967366238,
9.89919969516571e-05,
6.47728839110443e-05,
4.06093558895128e-05,
2.40706216722095e-05,
1.32278775917788e-05,
6.56367236603606e-06,
2.89543259446484e-06,
1.31177412220194e-06,
1.11918172661823e-06,
1.79836347095379e-06,
2.96751602061125e-06,
4.35308361161955e-06,
5.76499805329840e-06,
7.07765842414418e-06,
8.21354735367137e-06,
9.13114216166253e-06,
9.81469684998417e-06,
1.02666161123944e-05,
1.05013207681944e-05,
1.05407977677218e-05,
1.04111806786600e-05,
1.01400550619744e-05,
9.75498447394318e-06,
9.28203776880245e-06,
8.74513113628124e-06,
8.16560816526371e-06,
7.56192611136283e-06,
6.94975966207589e-06,
6.34202745021571e-06,
5.74905425138788e-06,
5.17892680624413e-06,
4.63757174402013e-06,
4.12921363490272e-06 };

  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[i] = temp_pp0[i];
      phiSq_pn0[i] = temp_pn0[i];
      phiSq_pn1[i] = temp_pn1[i];
    }
}


void Nuclear_Info::fill_arrays_chiral()
{
  double temp_pp0[100]={
503930.368696495,
317812.645921961,
81194.2990828164,
26342.2271162160,
10217.2745351813,
4463.56867060462,
2114.28516487998,
1058.82206862434,
550.633520907336,
293.283015533876,
158.153289161845,
85.4305981522729,
45.7238753124801,
23.9425897890113,
12.0641997554354,
5.70692890002787,
2.42915266130507,
0.852139388488781,
0.192368982636703,
0.00516991979812744,
0.0409256425824494,
0.162897538397258,
0.299575913689408,
0.416775763909872,
0.501207927970926,
0.550810474506498,
0.569102088852144,
0.561939449282667,
0.535715370645833,
0.496416905934596,
0.449188687661892,
0.398187399448207,
0.346597413311925,
0.296727704053509,
0.250143456852016,
0.207805779157984,
0.170204047007843,
0.137473493290365,
0.109495417862332,
0.0859794082410562,
0.0665286529342474,
0.0506903950435409,
0.0379935395042883,
0.0279755150874296,
0.0202005251250126,
0.0142709615194113,
0.00983352302026072,
0.00658142849063860,
0.00425381882565410,
0.00263321804127829,
0.00154178435676659,
0.000836895645053608,
0.000406466945880216,
0.000164305855940588,
4.57107776281637e-05,
3.43905976442862e-06,
4.12278794649675e-06,
2.51635486383012e-05,
5.21047910009569e-05,
7.64604671171241e-05,
9.39621644262655e-05,
0.000103179188959853,
0.000104462666364760,
9.91625979627999e-05,
8.90701024308095e-05,
7.60409167494763e-05,
6.17595449489965e-05,
4.76096861479643e-05,
3.46221337450286e-05,
2.34755146419335e-05,
1.45307103199138e-05,
7.88438064403275e-06,
3.43019788728449e-06,
9.19972101407957e-07,
1.96320453014268e-08,
3.56858024581210e-07,
1.55897152805634e-06,
3.28089650225942e-06,
5.22374207808405e-06,
7.14521096423176e-06,
8.86323384078690e-06,
1.02543147923701e-05,
1.12481604825284e-05,
1.18198990316728e-05,
1.19810052074559e-05,
1.17699689491139e-05,
1.12433863410400e-05,
1.04679253229323e-05,
9.51356204631237e-06,
8.44821860424595e-06,
7.33376752115819e-06,
6.22338677352309e-06,
5.16011461731331e-06,
4.17636806400599e-06,
3.29425783505339e-06,
2.52649760812261e-06,
1.87767748059471e-06,
1.34575052687488e-06,
9.23597826686550e-07,
6.00534757477642e-07
  };

  double temp_pn0[100]={
13608336.7628594,
620091.250856635,
111601.327671214,
31929.4599025633,
11626.3126988630,
4894.39350778088,
2263.92620734091,
1115.69668631523,
573.787872663755,
303.279258800782,
162.725310274508,
87.6583605305736,
46.8884310286590,
24.5960889806481,
12.4529922472042,
5.94568973854930,
2.57465764165572,
0.935118162346288,
0.231403645584911,
0.0127978086009854,
0.0263354305590984,
0.133144724286160,
0.260166313195804,
0.371984875049462,
0.454293299561301,
0.504179920408632,
0.524459935782691,
0.520418195987062,
0.497993613885044,
0.462822558442602,
0.419787582678229,
0.372857565167033,
0.325089534168373,
0.278713249880730,
0.235252032453815,
0.195653068921697,
0.160411989061703,
0.129684190325120,
0.103379966697558,
0.0812428125957117,
0.0629118312379605,
0.0479700297506312,
0.0359804272858151,
0.0265119700823780,
0.0191572159821438,
0.0135434868328238,
0.00933896862472178,
0.00625505540906111,
0.00404598314078754,
0.00250659466932110,
0.00146892205550713,
0.000798107827462472,
0.000388051934407082,
0.000157073926783409,
4.37881115347829e-05,
3.31733607938057e-06,
3.92031494382070e-06,
2.40645699288556e-05,
4.99461209291920e-05,
7.34367994047411e-05,
9.04249983755842e-05,
9.95076949500978e-05,
0.000100987862745088,
9.61299494977655e-05,
8.66284471955658e-05,
7.42479864208844e-05,
6.05970783631772e-05,
4.70029539357947e-05,
3.44600620451806e-05,
2.36291615584775e-05,
1.48687934198279e-05,
8.28513914055661e-06,
3.78956096158888e-06,
1.15633981345483e-06,
7.57110127819470e-08,
1.99170769147556e-07,
1.17566275836335e-06,
2.67841010135221e-06,
4.42290554632099e-06,
6.17714113195389e-06,
7.76539756250524e-06,
9.06700484680818e-06,
1.00115128737110e-05,
1.05715294275822e-05,
1.07543046701721e-05,
1.05930024576793e-05,
1.01383264075002e-05,
9.45096043650481e-06,
8.59516631310507e-06,
7.63367852112432e-06,
6.62389698061891e-06,
5.61533742191522e-06,
4.64820437136074e-06,
3.75288949057527e-06,
2.95022075111808e-06,
2.25227102179927e-06,
1.66352566181774e-06,
1.18225724997134e-06,
8.01974927730441e-07,
5.12829720284152e-07
  };

  double temp_pn1[100]={
2854560460.96972,
141323.717787954,
25009.3237360282,
7060.87092808851,
2556.81519903387,
1081.50618406432,
509.044752022738,
259.274413087024,
140.525864873376,
80.2358374102270,
47.9778025716885,
29.9462659853159,
19.4742895530370,
13.1751452834294,
9.25552780523634,
6.73276706470607,
5.05228863202439,
3.89308962181635,
3.06522875357505,
2.45398448413397,
1.98868365110418,
1.62494584093455,
1.33431072372797,
1.09809178636161,
0.903695514353953,
0.742375150751383,
0.607848220902178,
0.495439910234997,
0.401546656360873,
0.323297601798812,
0.258339235308859,
0.204698217965167,
0.160690534744196,
0.124863015848511,
0.0959528058027444,
0.0728590760913954,
0.0546222194866815,
0.0404076692550724,
0.0294928312295611,
0.0212555892482672,
0.0151642214399965,
0.0107680382577287,
0.00768878386247583,
0.00561243430248887,
0.00428152184905341,
0.00348801369735916,
0.00306657345183641,
0.00288836135723188,
0.00285539905586762,
0.00289539101327574,
0.00295708145992169,
0.00300618036074806,
0.00302175496234754,
0.00299309261510016,
0.00291708139358990,
0.00279599023284636,
0.00263562309662189,
0.00244388436154232,
0.00222961996284023,
0.00200176503962614,
0.00176871521127252,
0.00153792572999051,
0.00131564860073318,
0.00110682540080004,
0.000915078850830872,
0.000742781401747289,
0.000591174901257746,
0.000460530198804458,
0.000350316017131373,
0.000259378732488514,
0.000186113811262061,
0.000128622101583872,
8.48472643136566e-05,
5.26946353298040e-05,
3.01251481646567e-05,
1.52264209836987e-05,
6.26434276450128e-06,
1.71448033132140e-06,
2.77459830946586e-07,
8.80904903148086e-07,
2.67026396428592e-06,
4.99212090936746e-06,
7.37181441586429e-06,
9.48825061864297e-06,
1.11473113732220e-05,
1.22560491181754e-05,
1.27978857973236e-05,
1.28109064165820e-05,
1.23686845649394e-05,
1.15643227578629e-05,
1.04981090113555e-05,
9.26781360216352e-06,
7.96217614320125e-06,
6.65656822999872e-06,
5.41104673156757e-06,
4.26969648752344e-06,
3.26145791516896e-06,
2.40172225798250e-06,
1.69442147990566e-06,
1.13443105151659e-06
  };

  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[i] = temp_pp0[i];
      phiSq_pn0[i] = temp_pn0[i];
      phiSq_pn1[i] = temp_pn1[i];
    }
}


void Nuclear_Info::fill_arrays_chiral_n3lo()
{
  double temp_pp0[100]={
5554747.05824,
1571759.95008,
263239.844666,
45309.9688634,
12007.3480061,
4484.50966427,
1964.60077544,
929.823771006,
461.484841146,
234.760471422,
119.557836813,
59.9946952733,
29.086245682,
13.0830591515,
5.09860110795,
1.48711512766,
0.177944518449,
0.0385463160423,
0.44744347885,
1.03596695039,
1.6001599936,
2.034490268,
2.27946047573,
2.32198121127,
2.19735420468,
1.96079294838,
1.66386176718,
1.349557375,
1.04865970727,
0.777344096042,
0.54284541273,
0.351882885509,
0.209858869023,
0.115194701636,
0.0585342349351,
0.0276297518303,
0.0120803313581,
0.00486465773275,
0.00180713126043,
0.000631623314115,
0.000216363906923,
7.61359155349e-05,
2.81438084839e-05,
1.06849213989e-05,
3.9518898538e-06,
1.34971603446e-06,
4.10039510927e-07,
1.08554803195e-07,
2.48238592631e-08,
4.89416839249e-09,
8.33279489006e-10,
1.22913484412e-10,
1.57707504203e-11,
1.76714907627e-12,
1.73621048138e-13,
1.50131514941e-14,
1.14456935143e-15,
7.87958413981e-17,
2.55306871617e-18,
-8.77228159694e-19,
-1.16916156788e-18,
7.19060543055e-19,
4.7102027475e-19,
-1.22945216304e-19,
-2.20615980477e-20,
-3.55325795812e-19,
2.13674629059e-19,
1.67806178893e-20,
8.61362213598e-20,
2.49087420747e-21,
-8.62869478477e-20,
-7.68379836497e-20,
2.05321208439e-20,
1.41968485675e-21,
3.22705410637e-21,
-2.25320233494e-21,
2.65389680362e-21,
-8.76752181312e-22,
1.39763119167e-21,
2.86951501108e-23,
-4.11546775806e-22,
7.09453712693e-22,
2.34038570875e-22,
-6.25062745404e-23,
-1.03787086393e-22,
3.68653190533e-23,
6.78673777266e-23,
-1.95500187805e-23,
-1.85702966091e-24,
-4.95358636063e-24,
3.6531340888e-24,
1.4971106737e-24,
-2.2739073946e-24,
1.71471212442e-24,
2.37219693067e-25,
5.03291609111e-25,
-4.08294256856e-27,
2.8094624051e-25,
3.65242012122e-25,
8.07338667138e-26
  };

  double temp_pn0[100]={
3418811.00639,
1083955.11625,
218130.220233,
44642.9258224,
12718.3991313,
4786.49429673,
2091.81433622,
992.367661661,
494.061544089,
252.546447748,
130.214208212,
66.5399313852,
33.0055310344,
15.4297055628,
6.46389526189,
2.17778397297,
0.422312363427,
0.000108225274721,
0.231533539176,
0.726003651835,
1.25600788102,
1.69045246117,
1.96277476265,
2.05498641923,
1.98552976221,
1.79724774589,
1.54177889795,
1.26332806791,
0.990041090634,
0.737125796259,
0.515307109764,
0.334288534701,
0.199693563481,
0.109604193647,
0.0553718229971,
0.0258232360357,
0.0111334935509,
0.00443149360532,
0.00162173349741,
0.000542857115467,
0.000165430549179,
4.57602473719e-05,
1.14802093796e-05,
2.61557654557e-06,
5.42618783801e-07,
1.02834439413e-07,
1.78584187329e-08,
2.84749595506e-09,
4.17154185223e-10,
5.61027067735e-11,
6.91275323328e-12,
7.78266381501e-13,
7.98070608052e-14,
7.42913851011e-15,
6.31578915502e-16,
4.88801223741e-17,
5.11312217073e-18,
1.41943530929e-18,
-4.23186246203e-19,
4.75789449678e-19,
-2.04590711765e-19,
2.86412242607e-19,
2.67907491093e-19,
-1.11346284636e-19,
-2.6152928174e-19,
-4.40607156092e-20,
1.23138027551e-19,
5.4400337729e-20,
-7.47226519257e-21,
8.96467538823e-21,
-1.63024534263e-20,
-1.14109797051e-20,
1.14149518996e-21,
-2.93482748705e-21,
8.73826030531e-22,
6.20683753617e-22,
6.74989326136e-22,
1.77182570149e-22,
-1.88622490128e-23,
1.18036594981e-22,
-1.8102992419e-22,
1.04468713696e-22,
-1.09859549004e-22,
-5.39764478369e-23,
1.51681081863e-23,
-4.70103537069e-24,
3.16697387908e-23,
4.0686620163e-24,
-1.51181720277e-24,
-2.13352237781e-24,
-2.43126672369e-24,
-1.73931045274e-24,
-1.64284287356e-24,
-1.32495382732e-25,
2.29297761164e-25,
2.74609150565e-25,
-3.97548568394e-26,
7.97730134243e-26,
1.10080856979e-25,
1.8880407616e-26
  };

  double temp_pn1[100]={
42068.2689351,
23968.2042593,
11043.4021168,
4810.847213,
2177.29807512,
1051.05185441,
538.687509421,
290.017069217,
162.43575352,
93.9724129707,
55.8795838934,
34.0711631745,
21.3158565298,
13.7430031365,
9.20232147458,
6.4634283167,
4.79929070178,
3.76742265599,
3.09483370329,
2.61459645948,
2.22658017517,
1.87911440465,
1.55542936609,
1.26007062058,
1.00415810001,
0.793598465681,
0.623846156505,
0.48321625559,
0.361770163292,
0.256811360327,
0.170661421391,
0.105449312664,
0.06038148704,
0.0319884230423,
0.0156706292071,
0.0070964683199,
0.00296771207993,
0.0011432728139,
0.000404169718882,
0.000130506397091,
3.83136868471e-05,
1.01900179664e-05,
2.44944107366e-06,
5.31856327236e-07,
1.04529927994e-07,
1.87038281839e-08,
3.08199313528e-09,
4.75845769549e-10,
7.01531593004e-11,
9.9886778599e-12,
1.36821388441e-12,
1.7715361258e-13,
2.12100502494e-14,
2.30552398733e-15,
2.24810096716e-16,
1.95382082936e-17,
1.50868969851e-18,
1.03120033626e-19,
6.55202823601e-21,
4.88018628028e-22,
-6.75516602203e-23,
4.51217560094e-23,
1.38855334415e-22,
1.52199252698e-23,
-1.64957641079e-23,
-2.51348982692e-23,
2.10641813967e-23,
-1.78219379726e-24,
-5.58898563479e-24,
2.08523120244e-24,
-2.50822519404e-24,
-9.17022005454e-25,
1.4337778247e-25,
-1.58660626523e-25,
5.46057995462e-26,
-1.56272280872e-25,
1.41585239175e-25,
4.27462509263e-26,
-3.7764110536e-27,
-2.22010852889e-26,
6.08180663525e-27,
1.0588846009e-27,
-1.80291847404e-28,
5.18065557187e-28,
-1.09676427655e-27,
3.10972876219e-28,
2.36066138639e-29,
-3.98571231227e-29,
5.32421312229e-28,
-1.24812889232e-28,
-6.46240107126e-29,
-1.34900439466e-30,
-4.63544505389e-29,
1.62487115792e-29,
1.86260143418e-31,
8.21744990773e-30,
2.13769776336e-30,
2.69163986758e-30,
2.05556949033e-31,
2.82569539767e-31
  };

  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[i] = temp_pp0[i];
      phiSq_pn0[i] = temp_pn0[i];
      phiSq_pn1[i] = temp_pn1[i];
    }
}



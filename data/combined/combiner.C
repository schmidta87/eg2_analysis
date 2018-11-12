void combiner()
{
  TFile * f1p[4];

  //f1p[0] = new TFile("../SRC_e2p_C_GoodRuns_coulomb.root");
  //f1p[1] = new TFile("../SRC_e2p_Al_GoodRuns_coulomb.root");
  //f1p[2] = new TFile("../SRC_e2p_Fe_GoodRuns_coulomb.root");
  //f1p[3] = new TFile("../SRC_e2p_Pb_GoodRuns_coulomb.root");
  //TFile * outfile = new TFile("SRC_e2p_all.root","RECREATE");

  f1p[0] = new TFile("../SRC_e1p_C_GoodRuns_coulomb.root");
  f1p[1] = new TFile("../SRC_e1p_Al_GoodRuns_coulomb.root");
  f1p[2] = new TFile("../SRC_e1p_Fe_GoodRuns_coulomb.root");
  f1p[3] = new TFile("../SRC_e1p_Pb_GoodRuns_coulomb.root");
  TFile * outfile = new TFile("SRC_e1p_all.root","RECREATE");


  TTree * outtree = new TTree("T","Combined tree");

  const int maxP=2;
  Int_t nmb, target_type, ID[maxP];
  Float_t Q2, Xb, Pe[3], Pe_size, Ee, phi_e, theta_e, EC_in, EC_out, EC_tot, electron_z, Nu, q[3], q_size, theta_q, phi_q;
  Float_t Pp[maxP][3], Pp_size[maxP], pq_angle[maxP], Ep[maxP], Rp[maxP][3], phi_p[maxP], theta_p[maxP];
  Float_t TimeCorr4[maxP], SC_Edep[maxP], theta_Pmiss[maxP], phi_Pmiss[maxP], Pmiss[maxP][3], Pmiss_size[maxP];
  Float_t Pmiss_q_angle[maxP], pq_angle[maxP];

  outtree->Branch("nmb",&nmb,"nmb/I");
  outtree->Branch("target_type",&target_type,"target_type/I");
  outtree->Branch("ID",ID,"ID[nmb]/I");
  outtree->Branch("Q2",&Q2,"Q2/F");
  outtree->Branch("Xb",&Xb,"Xb/F");
  outtree->Branch("Pe",Pe,"Pe[3]/F");
  outtree->Branch("Pe_size",&Pe_size,"Pe_size/F");
  outtree->Branch("Ee",&Ee,"Ee/F");
  outtree->Branch("phi_e",&phi_e,"phi_e/F");
  outtree->Branch("theta_e",&theta_e,"theta_e/F");
  outtree->Branch("EC_in",&EC_in,"EC_in/F");
  outtree->Branch("EC_out",&EC_out,"EC_out/F");
  outtree->Branch("EC_tot",&EC_tot,"EC_tot/F");
  outtree->Branch("electron_z",&electron_z,"electron_z/F");
  outtree->Branch("Nu",&Nu,"Nu/F");
  outtree->Branch("q",q,"q[3]/F");
  outtree->Branch("q_size",&q_size,"q_size/F");
  outtree->Branch("theta_q",&theta_q,"theta_q/F");
  outtree->Branch("phi_q",&phi_q,"phi_q/F");
  outtree->Branch("Pp",Pp,"Pp[nmb][3]/F");
  outtree->Branch("Pmiss",Pmiss,"Pmiss[nmb][3]/F");
  outtree->Branch("Rp",Rp,"Rp[nmb][3]/F");
  outtree->Branch("Pp_size",Pp_size,"Pp_size[nmb]/F");
  outtree->Branch("pq_angle",pq_angle,"pq_angle[nmb]/F");
  outtree->Branch("Ep",Ep,"Ep[nmb]/F");
  outtree->Branch("phi_p",phi_p,"phi_p[nmb]/F");
  outtree->Branch("theta_p",theta_p,"theta_p[nmb]/F");
  outtree->Branch("TimeCorr4",TimeCorr4,"TimeCorr4[nmb]/F");
  outtree->Branch("SC_Edep",SC_Edep,"SC_Edep[nmb]/F");
  outtree->Branch("theta_Pmiss",theta_Pmiss,"theta_Pmiss[nmb]/F");
  outtree->Branch("phi_Pmiss",phi_Pmiss,"phi_Pmiss[nmb]/F");
  outtree->Branch("Pmiss_size",Pmiss_size,"Pmiss_size[nmb]/F");
  outtree->Branch("Pmiss_q_angle",Pmiss_q_angle,"Pmiss_q_angle[nmb]/F");
  outtree->Branch("pq_angle",pq_angle,"pq_angle[nmb]/F");

  for (int f=0 ; f<4 ; f++)
    {
      TTree * intree = (TTree*)f1p[f]->Get("T");
      // Set the branch addresses
      intree->SetBranchAddress("nmb",&nmb);
      intree->SetBranchAddress("target_type",&target_type);
      intree->SetBranchAddress("ID",ID);
      intree->SetBranchAddress("Q2",&Q2);
      intree->SetBranchAddress("Xb",&Xb);
      intree->SetBranchAddress("Pe",Pe);
      intree->SetBranchAddress("Pe_size",&Pe_size);
      intree->SetBranchAddress("Ee",&Ee);
      intree->SetBranchAddress("phi_e",&phi_e);
      intree->SetBranchAddress("theta_e",&theta_e);
      intree->SetBranchAddress("EC_in",&EC_in);
      intree->SetBranchAddress("EC_out",&EC_out);
      intree->SetBranchAddress("EC_tot",&EC_tot);
      intree->SetBranchAddress("electron_z",&electron_z);
      intree->SetBranchAddress("Nu",&Nu);
      intree->SetBranchAddress("q",q);
      intree->SetBranchAddress("q_size",&q_size);
      intree->SetBranchAddress("theta_q",&theta_q);
      intree->SetBranchAddress("phi_q",&phi_q);
      intree->SetBranchAddress("Pp",Pp);
      intree->SetBranchAddress("Pmiss",Pmiss);
      intree->SetBranchAddress("Rp",Rp);
      intree->SetBranchAddress("Pp_size",Pp_size);
      intree->SetBranchAddress("pq_angle",pq_angle);
      intree->SetBranchAddress("Ep",Ep);
      intree->SetBranchAddress("phi_p",phi_p);
      intree->SetBranchAddress("theta_p",theta_p);
      intree->SetBranchAddress("TimeCorr4",TimeCorr4);
      intree->SetBranchAddress("SC_Edep",SC_Edep);
      intree->SetBranchAddress("theta_Pmiss",theta_Pmiss);
      intree->SetBranchAddress("phi_Pmiss",phi_Pmiss);
      intree->SetBranchAddress("Pmiss_size",Pmiss_size);
      intree->SetBranchAddress("Pmiss_q_angle",Pmiss_q_angle);
      intree->SetBranchAddress("pq_angle",pq_angle);

      // Loop over events and fill them to the outtree
      for (int e=0 ; e<intree->GetEntries() ; e++)
	{
	  intree->GetEvent(e);
	  outtree->Fill();
	}
    }
  outtree->Write();
  outfile->Close();  
}

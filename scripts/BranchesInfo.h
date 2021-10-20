
double protonmass = 0.9382720;

double Jpsi_mass_Hepmc, Jpsi_theta_Hepmc, Jpsi_phi_Hepmc, Jpsi_eta_Hepmc, Jpsi_p_Hepmc, Jpsi_pT_Hepmc,Jpsi_E_Hepmc;
double ep_theta_Hepmc,  ep_phi_Hepmc,     ep_eta_Hepmc,   ep_p_Hepmc,     ep_pT_Hepmc,  ep_E_Hepmc;
double em_theta_Hepmc,  em_phi_Hepmc,     em_eta_Hepmc,   em_p_Hepmc,     em_pT_Hepmc,  em_E_Hepmc;
double e_theta_Hepmc,   e_phi_Hepmc,      e_eta_Hepmc,    e_p_Hepmc,      e_pT_Hepmc,   e_E_Hepmc;
double p_theta_Hepmc,   p_phi_Hepmc,      p_eta_Hepmc,    p_p_Hepmc,      p_pT_Hepmc,   p_E_Hepmc;
double Q2Hepmc, tHepmc;

double Jpsi_mass_MC, Jpsi_theta_MC, Jpsi_phi_MC, Jpsi_eta_MC, Jpsi_p_MC, Jpsi_pT_MC, Jpsi_E_MC;
double ep_theta_MC, ep_phi_MC, ep_pT_MC, ep_eta_MC, ep_p_MC, ep_E_MC;
double em_theta_MC, em_phi_MC, em_pT_MC, em_eta_MC, em_p_MC, em_E_MC;
double p_theta_MC, p_phi_MC, p_pT_MC, p_eta_MC, p_p_MC, p_E_MC;
double e_theta_MC, e_phi_MC, e_pT_MC, e_eta_MC, e_p_MC, e_E_MC;

void AddBranchesHepmc(TTree* inputTree){

	inputTree->Branch("Jpsi_mass_Hepmc",&Jpsi_mass_Hepmc,"Jpsi_mass_Hepmc/D");
	inputTree->Branch("Jpsi_theta_Hepmc",&Jpsi_theta_Hepmc,"Jpsi_theta_Hepmc/D");
	inputTree->Branch("Jpsi_eta_Hepmc",&Jpsi_eta_Hepmc,"Jpsi_eta_Hepmc/D");
	inputTree->Branch("Jpsi_p_Hepmc",&Jpsi_p_Hepmc,"Jpsi_p_Hepmc/D");
	inputTree->Branch("Jpsi_pT_Hepmc",&Jpsi_pT_Hepmc,"Jpsi_pT_Hepmc/D");
	inputTree->Branch("Jpsi_E_Hepmc",&Jpsi_E_Hepmc,"Jpsi_E_Hepmc/D");
	inputTree->Branch("ep_theta_Hepmc",&ep_theta_Hepmc,"ep_theta_Hepmc/D");
	inputTree->Branch("ep_phi_Hepmc",&ep_phi_Hepmc,"ep_phi_Hepmc/D");
	inputTree->Branch("ep_eta_Hepmc",&ep_eta_Hepmc,"ep_eta_Hepmc/D");
	inputTree->Branch("ep_p_Hepmc",&ep_p_Hepmc,"ep_p_Hepmc/D");
	inputTree->Branch("ep_pT_Hepmc",&ep_pT_Hepmc,"ep_pT_Hepmc/D");
	inputTree->Branch("ep_E_Hepmc",&ep_E_Hepmc,"ep_E_Hepmc/D");
	inputTree->Branch("em_theta_Hepmc",&em_theta_Hepmc,"em_theta_Hepmc/D");
	inputTree->Branch("em_phi_Hepmc",&em_phi_Hepmc,"em_phi_Hepmc/D");
	inputTree->Branch("em_eta_Hepmc",&em_eta_Hepmc,"em_eta_Hepmc/D");
	inputTree->Branch("em_p_Hepmc",&em_p_Hepmc,"em_p_Hepmc/D");
	inputTree->Branch("em_pT_Hepmc",&em_pT_Hepmc,"em_pT_Hepmc/D");
	inputTree->Branch("em_E_Hepmc",&em_E_Hepmc,"em_E_Hepmc/D");
	inputTree->Branch("e_theta_Hepmc",&e_theta_Hepmc,"e_theta_Hepmc/D");
	inputTree->Branch("e_phi_Hepmc",&e_phi_Hepmc,"e_phi_Hepmc/D");
	inputTree->Branch("e_eta_Hepmc",&e_eta_Hepmc,"e_eta_Hepmc/D");
	inputTree->Branch("e_p_Hepmc",&e_p_Hepmc,"e_p_Hepmc/D");
	inputTree->Branch("e_pT_Hepmc",&e_pT_Hepmc,"e_pT_Hepmc/D");
	inputTree->Branch("e_E_Hepmc",&e_E_Hepmc,"e_E_Hepmc/D");
	inputTree->Branch("p_theta_Hepmc",&p_theta_Hepmc,"p_theta_Hepmc/D");
	inputTree->Branch("p_phi_Hepmc",&p_phi_Hepmc,"p_phi_Hepmc/D");
	inputTree->Branch("p_eta_Hepmc",&p_eta_Hepmc,"p_eta_Hepmc/D");
	inputTree->Branch("p_p_Hepmc",&p_p_Hepmc,"p_p_Hepmc/D");
	inputTree->Branch("p_pT_Hepmc",&p_pT_Hepmc,"p_pT_Hepmc/D");
	inputTree->Branch("p_E_Hepmc",&p_E_Hepmc,"p_E_Hepmc/D");
	inputTree->Branch("tHepmc",&tHepmc,"tHepmc/D");

}

void AddBranchesMC(TTree* inputTree){

	inputTree->Branch("Jpsi_mass_MC",&Jpsi_mass_MC,"Jpsi_mass_MC/D");
	inputTree->Branch("Jpsi_theta_MC",&Jpsi_theta_MC,"Jpsi_theta_MC/D");
	inputTree->Branch("Jpsi_phi_MC",&Jpsi_phi_MC,"Jpsi_phi_MC/D");
	inputTree->Branch("Jpsi_eta_MC",&Jpsi_eta_MC,"Jpsi_eta_MC/D");
	inputTree->Branch("Jpsi_p_MC",&Jpsi_p_MC,"Jpsi_p_MC/D");
	inputTree->Branch("Jpsi_pT_MC",&Jpsi_pT_MC,"Jpsi_pT_MC/D");
	inputTree->Branch("Jpsi_E_MC",&Jpsi_E_MC,"Jpsi_E_MC/D");

	inputTree->Branch("ep_theta_MC",&ep_theta_MC,"ep_theta_MC/D");
	inputTree->Branch("ep_phi_MC",&ep_phi_MC,"ep_phi_MC/D");
	inputTree->Branch("ep_eta_MC",&ep_eta_MC,"ep_eta_MC/D");
	inputTree->Branch("ep_pT_MC",&ep_pT_MC,"ep_pT_MC/D");
	inputTree->Branch("ep_p_MC",&ep_p_MC,"ep_p_MC/D");
	inputTree->Branch("ep_E_MC",&ep_E_MC,"ep_E_MC/D");

	inputTree->Branch("em_theta_MC",&em_theta_MC,"em_theta_MC/D");
	inputTree->Branch("em_phi_MC",&em_phi_MC,"em_phi_MC/D");
	inputTree->Branch("em_eta_MC",&em_eta_MC,"em_eta_MC/D");
	inputTree->Branch("em_pT_MC",&em_pT_MC,"em_pT_MC/D");
	inputTree->Branch("em_p_MC",&em_p_MC,"em_p_MC/D");
	inputTree->Branch("em_E_MC",&em_E_MC,"em_E_MC/D");

	inputTree->Branch("p_theta_MC",&p_theta_MC,"p_theta_MC/D");
	inputTree->Branch("p_phi_MC",  &p_phi_MC,  "p_phi_MC/D");
	inputTree->Branch("p_eta_MC",  &p_eta_MC,  "p_eta_MC/D");
	inputTree->Branch("p_pT_MC",   &p_pT_MC,   "p_pT_MC/D");
	inputTree->Branch("p_p_MC",    &p_p_MC,    "p_p_MC/D");
	inputTree->Branch("p_E_MC",    &p_E_MC,    "p_E_MC/D");

	inputTree->Branch("e_theta_MC",&e_theta_MC,"e_theta_MC/D");
	inputTree->Branch("e_phi_MC",&e_phi_MC,"e_phi_MC/D");
	inputTree->Branch("e_eta_MC",&e_eta_MC,"e_eta_MC/D");
	inputTree->Branch("e_pT_MC",&e_pT_MC,"e_pT_MC/D");
	inputTree->Branch("e_p_MC",&e_p_MC,"e_p_MC/D");
	inputTree->Branch("e_E_MC",&e_E_MC,"e_E_MC/D");


}



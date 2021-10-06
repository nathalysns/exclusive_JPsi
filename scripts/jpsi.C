

#include "jpsi.h"
#include "BranchesInfo.h"
//#include "caloheader.h"
//#include "clusterizer.cxx"

void jpsi(
	TString inFile            	= "./../rootfiles/JPsi_all.root",
	TString inFileGeometry      = "./../rootfiles/geometry.root",
	bool do_reclus              = true,
    unsigned short primaryTrackSource = 0

){


	// load tree
    TChain *const tt_event = new TChain("event_tree");
    if (inFile.EndsWith(".root")) {                     // are we loading a single root tree?
        std::cout << "loading a single root file" << std::endl;
        tt_event->AddFile(inFile);
    }
    else {                                              // or are we loading a bunch?
        std::cout << "loading a list of files" << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
            tt_event->AddFile(filePath.c_str());
        }
        files.close();
    }
    if(!tt_event){ std::cout << "tree not found... returning!"<< std::endl; return;}

     // // load geometry tree
    /*tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("geometry_tree");
    if(!tt_geometry){ cout << "geometry tree not found... returning!"<< endl; return;}*/

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    std::cout << "Number of events in tree: " << nEntriesTree << std::endl;

    SetBranchAddressesTree(tt_event);
    //SetBranchAddressesGeometryTree(tt_geometry);
    //SetGeometryIndices();

    double elect_energy = 10;
    double proton_energy = 100;

    TFile *MyFile = new TFile("jpsi_test.root","RECREATE");
	TTree *EvTree = new TTree( "T", "T" );

	int nTracks;
	int n_mass;
	int trackID[20];
	

	double Jpsi_track3_M1, Jpsi_track3_M2, Jpsi_M;
	double ep_eta, em_eta, p_eta, Jpsi_eta, e_eta;
	double ep_theta, em_theta, p_theta, Jpsi_theta, e_theta;
	double ep_phi, em_phi, p_phi, Jpsi_phi, e_phi;
	double ep_p, em_p, p_p, Jpsi_p, e_p;
	double ep_pT, em_pT, p_pT, Jpsi_pT, e_pT;
	double ep_E, em_E, p_E, Jpsi_E, e_E, phiJpsidet;

	double Q2, xbj, Q2MC, xbjMC, phiJpsi;
	double ZDCsurrogate_hits_x[20], ZDCsurrogate_hits_y[20], ZDCsurrogate_hits_z[20];
    double ZDCsurrogate_hits_x2[20], ZDCsurrogate_hits_y2[20], ZDCsurrogate_hits_z2[20];
    double ZDCsurrogate_hits_t[20], ZDCsurrogate_hits_edep[20];
    double rp1_hits_x[20], rp1_hits_y[20], rp1_hits_z[20];
    double rp1_hits_x2[20], rp1_hits_y2[20], rp1_hits_z2[20];
    double rp1_hits_t[20], rp1_hits_edep[20];
    double rp2_hits_x[20], rp2_hits_y[20], rp2_hits_z[20];
    double rp2_hits_x2[20], rp2_hits_y2[20], rp2_hits_z2[20];
    double rp2_hits_t[20], rp2_hits_edep[20];
    double b0_hits_x[20], b0_hits_y[20], b0_hits_z[20];
    double b0_hits_x2[20], b0_hits_y2[20], b0_hits_z2[20];
    double b0_hits_t[20], b0_hits_edep[20];
    double offmom_hits_x[20], offmom_hits_y[20], offmom_hits_z[20];
    double offmom_hits_x2[20], offmom_hits_y2[20], offmom_hits_z2[20];
    double offmom_hits_t[20], offmom_hits_edep[20];

	int nHits; 
	int hits_layerID[20];
	int hits_trueID[20];

	EvTree->Branch("Jpsi_track3_M1",&Jpsi_track3_M1,"Jpsi_track3_M1/D");
	EvTree->Branch("Jpsi_track3_M2",&Jpsi_track3_M2,"Jpsi_track3_M2/D"); 	
	EvTree->Branch("Jpsi_M",&Jpsi_M,"Jpsi_M/D");

	EvTree->Branch("ep_eta",&ep_eta,"ep_eta/D");
	EvTree->Branch("ep_theta",&ep_theta,"ep_theta/D"); 	
	EvTree->Branch("ep_phi",&ep_phi,"ep_phi/D");
	EvTree->Branch("ep_p",&ep_p,"ep_p/D");
	EvTree->Branch("ep_pT",&ep_pT,"ep_pT/D");
	EvTree->Branch("ep_E",&ep_E,"ep_E/D");
	EvTree->Branch("em_eta",&em_eta,"em_eta/D");
	EvTree->Branch("em_theta",&em_theta,"em_theta/D"); 	
	EvTree->Branch("em_phi",&em_phi,"em_phi/D");
	EvTree->Branch("em_p",&em_p,"em_p/D");
	EvTree->Branch("em_pT",&em_pT,"em_pT/D");
	EvTree->Branch("em_E",&em_E,"em_E/D");
	EvTree->Branch("p_eta",&p_eta,"p_eta/D");
	EvTree->Branch("p_theta",&p_theta,"p_theta/D"); 	
	EvTree->Branch("p_phi",&p_phi,"p_phi/D");
	EvTree->Branch("p_p",&p_p,"p_p/D");
	EvTree->Branch("p_pT",&p_pT,"p_pT/D");
	EvTree->Branch("p_E",&p_E,"p_E/D");
	EvTree->Branch("Jpsi_eta",&Jpsi_eta,"Jpsi_eta/D");
	EvTree->Branch("Jpsi_theta",&Jpsi_theta,"Jpsi_theta/D"); 	
	EvTree->Branch("Jpsi_phi",&Jpsi_phi,"Jpsi_phi/D");
	EvTree->Branch("Jpsi_p",&Jpsi_p,"Jpsi_p/D");
	EvTree->Branch("Jpsi_pT",&Jpsi_pT,"Jpsi_pT/D");
	EvTree->Branch("Jpsi_E",&Jpsi_E,"Jpsi_E/D");
	EvTree->Branch("e_eta",&e_eta,"e_eta/D");
	EvTree->Branch("e_theta",&e_theta,"e_theta/D"); 	
	EvTree->Branch("e_phi",&e_phi,"e_phi/D");
	EvTree->Branch("e_p",&e_p,"e_p/D");
	EvTree->Branch("e_pT",&e_pT,"e_pT/D");
	EvTree->Branch("e_E",&e_E,"e_E/D");
	
	EvTree->Branch("Q2",&Q2,"Q2/D");
	EvTree->Branch("xbj",&xbj,"xbj/D");
	EvTree->Branch("Q2MC",&Q2MC,"Q2MC/D");
	EvTree->Branch("xbjMC",&xbjMC,"xbjMC/D");
	EvTree->Branch("phiJpsi",&phiJpsi,"phiJpsi/D");
	EvTree->Branch("phiJpsidet",&phiJpsidet,"phiJpsidet/D");

	EvTree->Branch("nHits",&nHits,"nHits/I");
	EvTree->Branch("hits_layerID",hits_layerID,"hits_layerID[nHits]/I");
	EvTree->Branch("hits_trueID",hits_trueID,"hits_trueID[nHits]/I");
	EvTree->Branch("ZDCsurrogate_hits_x",    ZDCsurrogate_hits_x,    "ZDCsurrogate_hits_x[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_y",    ZDCsurrogate_hits_y,    "ZDCsurrogate_hits_y[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_z",    ZDCsurrogate_hits_z,    "ZDCsurrogate_hits_z[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_x2",   ZDCsurrogate_hits_x2,   "ZDCsurrogate_hits_x2[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_y2",   ZDCsurrogate_hits_y2,   "ZDCsurrogate_hits_y2[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_z2",   ZDCsurrogate_hits_z2,   "ZDCsurrogate_hits_z2[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_t",    ZDCsurrogate_hits_t,    "ZDCsurrogate_hits_t[nHits]/D");
	EvTree->Branch("ZDCsurrogate_hits_edep", ZDCsurrogate_hits_edep, "ZDCsurrogate_hits_edep[nHits]/D");

	EvTree->Branch("rp1_hits_x",  rp1_hits_x,  "rp1_hits_x[nHits]/D");
	EvTree->Branch("rp1_hits_y",  rp1_hits_y,  "rp1_hits_y[nHits]/D");
	EvTree->Branch("rp1_hits_z",  rp1_hits_z,  "rp1_hits_z[nHits]/D");
	EvTree->Branch("rp1_hits_x2", rp1_hits_x2, "rp1_hits_x2[nHits]/D");
	EvTree->Branch("rp1_hits_y2", rp1_hits_y2, "rp1_hits_y2[nHits]/D");
	EvTree->Branch("rp1_hits_z2", rp1_hits_z2, "rp1_hits_z2[nHits]/D");
	EvTree->Branch("rp1_hits_t",  rp1_hits_t,  "rp1_hits_t[nHits]/D");
	EvTree->Branch("rp1_hits_edep",  rp1_hits_edep,  "rp1_hits_edep[nHits]/D");

	EvTree->Branch("rp2_hits_x", rp2_hits_x,  "rp2_hits_x[nHits]/D");
	EvTree->Branch("rp2_hits_y", rp2_hits_y,  "rp2_hits_y[nHits]/D");
	EvTree->Branch("rp2_hits_z", rp2_hits_z,  "rp2_hits_z[nHits]/D");
	EvTree->Branch("rp2_hits_x2",rp2_hits_x2, "rp2_hits_x2[nHits]/D");
	EvTree->Branch("rp2_hits_y2",rp2_hits_y2, "rp2_hits_y2[nHits]/D");
	EvTree->Branch("rp2_hits_z2",rp2_hits_z2, "rp2_hits_z2[nHits]/D");
	EvTree->Branch("rp2_hits_t", rp2_hits_t,  "rp2_hits_t[nHits]/D");
	EvTree->Branch("rp2_hits_edep", rp2_hits_edep,  "rp2_hits_edep[nHits]/D");

    EvTree->Branch("b0_hits_x",     b0_hits_x,  "b0_hits_x[nHits]/D");
	EvTree->Branch("b0_hits_y",     b0_hits_y,  "b0_hits_y[nHits]/D");
	EvTree->Branch("b0_hits_z",     b0_hits_z,  "b0_hits_z[nHits]/D");
	EvTree->Branch("b0_hits_x2",    b0_hits_x2, "b0_hits_x2[nHits]/D");
	EvTree->Branch("b0_hits_y2",    b0_hits_y2,   "b0_hits_y2[nHits]/D");
	EvTree->Branch("b0_hits_z2",    b0_hits_z2,    "b0_hits_z2[nHits]/D");
	EvTree->Branch("b0_hits_t",     b0_hits_t,     "b0_hits_t[nHits]/D");
	EvTree->Branch("b0_hits_edep",  b0_hits_edep,  "b0_hits_edep[nHits]/D");


	EvTree->Branch("offmom_hits_x",     offmom_hits_x,     "offmom_hits_x[nHits]/D");
	EvTree->Branch("offmom_hits_y",     offmom_hits_y,     "offmom_hits_y[nHits]/D");
	EvTree->Branch("offmom_hits_z",     offmom_hits_z,     "offmom_hits_z[nHits]/D");
	EvTree->Branch("offmom_hits_x2",    offmom_hits_x2,    "offmom_hits_x2[nHits]/D");
	EvTree->Branch("offmom_hits_y2",    offmom_hits_y2,    "offmom_hits_y2[nHits]/D");
	EvTree->Branch("offmom_hits_z2",    offmom_hits_z2,    "offmom_hits_z2[nHits]/D");
	EvTree->Branch("offmom_hits_t",     offmom_hits_t,     "offmom_hits_t[nHits]/D");
	EvTree->Branch("offmom_hits_edep",  offmom_hits_edep,  "offmom_hits_edep[nHits]/D");

	if(HepmcEnabled){
    		AddBranchesHepmc(EvTree);
    		cout << "test" << endl;
	}

	AddBranchesMC(EvTree);

   	_nEventsTree=0;

    int track = 0;
    int MC_aboveQ2 = 0;

    TLorentzVector ProtonBeam(0,0,274.998,275);
    TLorentzVector ElectronBeam(0,0,-18,18);
    TLorentzRotation rotlor = TLorentzRotation().RotateY(12.5e-3).Boost(sin(12.5e-3),0,0);

    int notaccepted =0 ;
    int alltrack3 = 0 ;

    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {  

    	tt_event->GetEntry(i);
    	
    	if (_hepmcp_Q2 < 1.) continue;

    	nHits = _nHits;

    	for(int h=0; h<nHits; h++){
    		hits_trueID[h]  = _hits_trueID[h];
    		hits_layerID[h] = _hits_layerID[h];
    		ZDCsurrogate_hits_x[h] = ZDCsurrogate_hits_y[h] = ZDCsurrogate_hits_z[h] = -1000;
    		ZDCsurrogate_hits_x2[h] = ZDCsurrogate_hits_y2[h] = ZDCsurrogate_hits_z2[h] = -1000;
    		ZDCsurrogate_hits_t[h] = ZDCsurrogate_hits_edep[h] = -1000;
    		rp1_hits_x[h]     = rp1_hits_y[h] = rp1_hits_z[h] = -1000;
    		rp1_hits_x2[h]    = rp1_hits_y2[h] = rp1_hits_z2[h] = -1000;
    		rp1_hits_t[h]     = rp1_hits_edep[h] = -1000;
    		rp2_hits_x[h]     = rp2_hits_y[h] = rp2_hits_z[h] = -1000;
    		rp2_hits_x2[h]    = rp2_hits_y2[h] = rp2_hits_z2[h] = -1000;
    		rp2_hits_t[h]     = rp2_hits_edep[h] = -1000;
    		b0_hits_x[h]      = b0_hits_y[h] = b0_hits_z[h] = -1000;
    		b0_hits_x2[h]     = b0_hits_y2[h] = b0_hits_z2[h] = -1000;
    		b0_hits_t[h]      = b0_hits_edep[h] = -1000;
    		offmom_hits_x[h]  = offmom_hits_y[h] = offmom_hits_z[h] = -1000;
    		offmom_hits_x2[h] = offmom_hits_y2[h] = offmom_hits_z2[h] = -1000;
    		offmom_hits_t[h]  = offmom_hits_edep[h] = -1000;

    		if(hits_layerID[h] == 70){
    			ZDCsurrogate_hits_x[h]    = _hits_x[h];
    			ZDCsurrogate_hits_y[h]    = _hits_y[h];
    			ZDCsurrogate_hits_z[h]    = _hits_z[h];
    			ZDCsurrogate_hits_x2[h]   = _hits_x2[h];
    			ZDCsurrogate_hits_y2[h]   = _hits_y2[h];
    			ZDCsurrogate_hits_z2[h]   = _hits_z2[h];
    			ZDCsurrogate_hits_t[h]    = _hits_t[h];
    			ZDCsurrogate_hits_edep[h] = _hits_edep[h]; 
    		}
    		if(hits_layerID[h] == 71){
    			rp1_hits_x[h]    = _hits_x[h];
    			rp1_hits_y[h]    = _hits_y[h];
    			rp1_hits_z[h]    = _hits_z[h];
    			rp1_hits_x2[h]   = _hits_x2[h];
    			rp1_hits_y2[h]   = _hits_y2[h];
    			rp1_hits_z2[h]   = _hits_z2[h];
    			rp1_hits_t[h]    = _hits_t[h];
    			rp1_hits_edep[h] = _hits_edep[h]; 
    		}

    		if(hits_layerID[h] == 72){
    			b0_hits_x[h]    = _hits_x[h];
    			b0_hits_y[h]    = _hits_y[h];
    			b0_hits_z[h]    = _hits_z[h];
    			b0_hits_x2[h]   = _hits_x2[h];
    			b0_hits_y2[h]   = _hits_y2[h];
    			b0_hits_z2[h]   = _hits_z2[h];
    			b0_hits_t[h]    = _hits_t[h];
    			b0_hits_edep[h] = _hits_edep[h]; 
    		}

    		if(hits_layerID[h] == 73){
    			rp2_hits_x[h]    = _hits_x[h];
    			rp2_hits_y[h]    = _hits_y[h];
    			rp2_hits_z[h]    = _hits_z[h];
    			rp2_hits_x2[h]   = _hits_x2[h];
    			rp2_hits_y2[h]   = _hits_y2[h];
    			rp2_hits_z2[h]   = _hits_z2[h];
    			rp2_hits_t[h]    = _hits_t[h];
    			rp2_hits_edep[h] = _hits_edep[h]; 
    		}
    		if(hits_layerID[h] == 74){
    			offmom_hits_x[h]    = _hits_x[h];
    			offmom_hits_y[h]    = _hits_y[h];
    			offmom_hits_z[h]    = _hits_z[h];
    			offmom_hits_x2[h]   = _hits_x2[h];
    			offmom_hits_y2[h]   = _hits_y2[h];
    			offmom_hits_z2[h]   = _hits_z2[h];
    			offmom_hits_t[h]    = _hits_t[h];
    			offmom_hits_edep[h] = _hits_edep[h]; 
    		}


    	}

    
    	TLorentzVector JpsiHepmc;
    	TLorentzVector eeHepmc;
    	TLorentzVector emHepmc;
    	TLorentzVector epHepmc;
    	TLorentzVector pHepmc;

    	for(int j = 0; j < _nHepmcp; j++){
    		if (_hepmcp_BCID[j] == 10006){
    			JpsiHepmc.SetPxPyPzE(_hepmcp_px[j],_hepmcp_py[j],_hepmcp_pz[j],_hepmcp_E[j]);
	    		Jpsi_mass_Hepmc = JpsiHepmc.M(); Jpsi_theta_Hepmc = JpsiHepmc.Theta();
	    		Jpsi_phi_Hepmc = JpsiHepmc.Phi(); Jpsi_eta_Hepmc  = JpsiHepmc.Eta(); 
	    		Jpsi_p_Hepmc = JpsiHepmc.Vect().Mag(); Jpsi_pT_Hepmc = JpsiHepmc.Pt();
	    		Jpsi_E_Hepmc = JpsiHepmc.E();
    		}
    		else if(_hepmcp_BCID[j] == 10002){
    			eeHepmc.SetPxPyPzE(_hepmcp_px[j],_hepmcp_py[j],_hepmcp_pz[j],_hepmcp_E[j]);
	    		e_theta_Hepmc = eeHepmc.Theta();
	    		e_phi_Hepmc = eeHepmc.Phi(); e_eta_Hepmc  = eeHepmc.Eta(); 
	    		e_p_Hepmc = eeHepmc.Vect().Mag(); e_pT_Hepmc = eeHepmc.Pt();
	    		e_E_Hepmc = eeHepmc.E();
    		}else if(_hepmcp_BCID[j] == 10009){
    			emHepmc.SetPxPyPzE(_hepmcp_px[j],_hepmcp_py[j],_hepmcp_pz[j],_hepmcp_E[j]);
	    		em_theta_Hepmc = emHepmc.Theta();
	    		em_phi_Hepmc = emHepmc.Phi(); em_eta_Hepmc  = emHepmc.Eta(); 
	    		em_p_Hepmc = emHepmc.Vect().Mag(); em_pT_Hepmc = emHepmc.Pt();
	    		em_E_Hepmc = emHepmc.E();
    		}else if(_hepmcp_BCID[j] == 10008){
    			epHepmc.SetPxPyPzE(_hepmcp_px[j],_hepmcp_py[j],_hepmcp_pz[j],_hepmcp_E[j]);
	    		ep_theta_Hepmc = epHepmc.Theta();
	    		ep_phi_Hepmc = epHepmc.Phi(); ep_eta_Hepmc  = epHepmc.Eta(); 
	    		ep_p_Hepmc = epHepmc.Vect().Mag(); ep_pT_Hepmc = epHepmc.Pt();
	    		ep_E_Hepmc = epHepmc.E();
    		}else if(_hepmcp_BCID[j] == 10007){
    			pHepmc.SetPxPyPzE(_hepmcp_px[j],_hepmcp_py[j],_hepmcp_pz[j],_hepmcp_E[j]);
	    		p_theta_Hepmc = pHepmc.Theta();
	    		p_phi_Hepmc = pHepmc.Phi(); p_eta_Hepmc  = pHepmc.Eta(); 
	    		p_p_Hepmc = pHepmc.Vect().Mag(); p_pT_Hepmc = pHepmc.Pt();
	    		p_E_Hepmc = pHepmc.E();
    		}else continue;


    	}
    	tHepmc = (pHepmc - ProtonBeam)*(pHepmc - ProtonBeam);
    	
    	// load current event
        TLorentzVector epMC_det;
        TLorentzVector emMC_det;
        TLorentzVector pMC_det;
        TLorentzVector eMC_det;

        for(int imc=0; imc<_nMCPart; imc++){
        	 //===== J/Psi True
        	 if(_mcpart_PDG[imc] == 11 && _mcpart_BCID[imc]==10009)  {
        	 	epMC_det[0] = _mcpart_px[imc];
        	 	epMC_det[1]	= _mcpart_py[imc];
        	 	epMC_det[2]	= _mcpart_pz[imc];
        	 	epMC_det[3]	= _mcpart_E[imc];
        	 }
        	 else if(_mcpart_PDG[imc] == -11 && _mcpart_BCID[imc]==10008)  {
        	 	emMC_det[0] = _mcpart_px[imc];
        	 	emMC_det[1]	= _mcpart_py[imc];
        	 	emMC_det[2]	= _mcpart_pz[imc];
        	 	emMC_det[3]	= _mcpart_E[imc];
        	 }
        	 else if(_mcpart_PDG[imc] == 2212 && _mcpart_BCID[imc]==10007)  {
        	 	pMC_det[0] = _mcpart_px[imc];
        	 	pMC_det[1]	= _mcpart_py[imc];
        	 	pMC_det[2]	= _mcpart_pz[imc];
        	 	pMC_det[3]	= _mcpart_E[imc];
        	 }
        	 else if(_mcpart_PDG[imc] == 11 && _mcpart_BCID[imc]==10002)  {
        	 	eMC_det[0] = _mcpart_px[imc];
        	 	eMC_det[1]	= _mcpart_py[imc];
        	 	eMC_det[2]	= _mcpart_pz[imc];
        	 	eMC_det[3]	= _mcpart_E[imc];
        	 }
        	 else
        	 {

        	 	continue;
        	 }
        }

        TLorentzVector JpsiMC_det 	= 	epMC_det + emMC_det;
/*
        TLorentzVector pMC        	=  	rotlor*pMC_det;  	
       	TLorentzVector JpsiMC 		= 	rotlor*JpsiMC_det;
       	TLorentzVector eMC    		=	rotlor*eMC_det;
       	TLorentzVector epMC    		=	rotlor*epMC_det;
       	TLorentzVector emMC    		=	rotlor*emMC_det;
*/
        TLorentzVector pMC        	=  	pMC_det;  	
       	TLorentzVector JpsiMC 		= 	JpsiMC_det;
       	TLorentzVector eMC    		=	eMC_det;
       	TLorentzVector epMC    		=	epMC_det;
       	TLorentzVector emMC    		=	emMC_det;

        Jpsi_mass_MC  = JpsiMC.M();
      	Jpsi_theta_MC = JpsiMC.Theta()*180/TMath::Pi();
      	Jpsi_phi_MC   = JpsiMC.Phi()*180/TMath::Pi();
        Jpsi_eta_MC   = JpsiMC.Eta();
        Jpsi_p_MC 	  = JpsiMC.Vect().Mag();
        ep_pT_MC	  = JpsiMC.Perp();
        Jpsi_E_MC     = JpsiMC.E();
        //====== ep
        ep_phi_MC     = epMC.Phi()*180/TMath::Pi();
      	ep_theta_MC   = epMC.Theta()*180/TMath::Pi();
        ep_eta_MC 	  = epMC.Eta();
        ep_pT_MC	  = epMC.Perp();
        ep_p_MC       = epMC.Vect().Mag();
        ep_E_MC       = epMC.E();

        //====== em
       	em_phi_MC     = emMC.Phi()*180/TMath::Pi();
      	em_theta_MC   = emMC.Theta()*180/TMath::Pi();
      	em_eta_MC 	  = -TMath::Log(TMath::Tan(emMC.Theta()/2));
        em_pT_MC	  = emMC.Perp();
        em_p_MC       = emMC.Vect().Mag();
        em_E_MC       = emMC.E();

        //====== p
        p_phi_MC     =  pMC.Phi();
        p_theta_MC   =  pMC.Theta()*180/TMath::Pi();
        p_eta_MC 	 =  pMC.Eta();
        p_pT_MC		 =  pMC.Perp();
        p_p_MC       =  pMC.Vect().Mag();
        p_E_MC       =  pMC.E();

        //====== e
      	e_phi_MC     =  eMC.Phi()*180/TMath::Pi();
      	e_theta_MC   =  eMC.Theta()*180/TMath::Pi();
      	e_eta_MC  	 =  eMC.Eta();
      	e_pT_MC	  	 =  eMC.Perp();
      	e_p_MC       =  eMC.Vect().Mag();
        e_E_MC       =  eMC.E();


		TLorentzVector hypep;
		TLorentzVector hyp1;
		TLorentzVector hyp2;

		int eptrack_N[3];
		int othercharged[3];
		int noe = 0;
		int nep = 0;
		if(_nTracks != 9) continue;

		alltrack3++;
		
		if(_nTracks ==9){
			
	
			//======= Read from the tracking paths ==============//
			for(int l=0; l<_nTracks; l++){
				if(_track_source[l]!=0) continue;
				if( _track_charge[l] == short(1)) {
					eptrack_N[nep] = l;
					nep++; 

				}
				else if( _track_charge[l] == short(-1))  {
					othercharged[noe] = l;
					noe++;
				}
			}
		}

		if(nep>1) continue;
		xsif(nep>2) continue;
			
		hypep.SetXYZM(_track_px[eptrack_N[0]],_track_py[eptrack_N[0]],_track_pz[eptrack_N[0]], me);
		hyp1.SetXYZM(_track_px[othercharged[0]],_track_py[othercharged[0]],_track_pz[othercharged[0]], me);
		hyp2.SetXYZM(_track_px[othercharged[1]],_track_py[othercharged[1]],_track_pz[othercharged[1]], me);

		TLorentzVector JPsihyp1 = hypep + hyp1;
		TLorentzVector JPsihyp2 = hypep + hyp2;
		
		TLorentzVector JPsi;
		TLorentzVector Electron;
		TLorentzVector ep;
		TLorentzVector em;
	
		ep = hypep;
		Jpsi_track3_M1 = JPsihyp1.M();
		Jpsi_track3_M2 = JPsihyp2.M();

		if(Jpsi_track3_M1 > low_limjpsiM && Jpsi_track3_M1 < high_limjpsiM){
			JPsi = JPsihyp1;
			Electron.SetXYZM(_track_px[othercharged[1]],_track_py[othercharged[1]],_track_pz[othercharged[1]], me);
			em = hyp1;
		} else if(Jpsi_track3_M2 > low_limjpsiM && Jpsi_track3_M2 < high_limjpsiM){
			JPsi = JPsihyp2;
			Electron.SetXYZM(_track_px[othercharged[0]],_track_py[othercharged[0]],_track_pz[othercharged[0]], me);
			em = hyp2;
		}else{
			notaccepted++;
			continue;
		}



		Jpsi_M = JPsi.M();
		ep_eta = hypep.Eta(); em_eta = em.Eta(); Jpsi_eta = JPsi.Eta(); e_eta = Electron.Eta();
		ep_theta = hypep.Theta(); em_theta = em.Theta(); Jpsi_theta = JPsi.Theta(); e_theta = Electron.Theta();
		ep_phi = hypep.Phi(); em_phi = em.Phi(); Jpsi_phi = JPsi.Phi(); e_phi = Electron.Phi();
		ep_p = hypep.Vect().Mag(); em_p = em.Vect().Mag(); Jpsi_p = JPsi.Vect().Mag(); e_p = Electron.Vect().Mag();
		ep_pT = hypep.Perp(); em_pT = em.Perp(); Jpsi_pT = JPsi.Perp(); e_pT = Electron.Perp();
		ep_E = hypep.E(); em_E = em.E(); Jpsi_E = JPsi.E(); e_E = Electron.E();

		//==== Kinematics from tracking ==========

	   	
    	TLorentzVector eMC_HO = rotate_reco(eMC,25.e-3);
    	TLorentzVector Electron_HO = rotate_reco(Electron,25.e-3);

		Q2MC = -(ElectronBeam - eMC_HO)*(ElectronBeam - eMC_HO);
		Q2   = -(ElectronBeam - Electron_HO)*(ElectronBeam - Electron_HO);


		xbjMC = Q2MC/(2*ProtonBeam*(ElectronBeam - eMC_HO));
		xbj   = Q2/(2*ProtonBeam*(ElectronBeam - Electron_HO));

		TLorentzVector JPsi_HO = rotate_reco(JPsi,25.e-3);

		phiJpsi = JPsi_HO.Vect().Angle(ProtonBeam.Vect());


		TLorentzVector q = (ElectronBeam - Electron_HO);
		q.Boost( -ProtonBeam.BoostVector() );
		Electron_HO.Boost( -ProtonBeam.BoostVector() );
		JPsi_HO.Boost( -ProtonBeam.BoostVector() );

		TVector3 v1      =  q.Vect().Cross(Electron_HO.Vect());
		TVector3 v2      =  q.Vect().Cross(JPsi_HO.Vect());

		phiJpsi          =  v1.Angle(v2);

		if(q.Vect().Dot(v1.Cross(v2))<0) phiJpsi=2.*TMath::Pi()-phiJpsi;
		
		TLorentzVector test = rotate_reco(ProtonBeam,25.e-3);


		phiJpsidet = JPsi.Vect().Angle(ProtonBeam.Vect());

		track++;
	


	  	EvTree->Fill();
    	//if (i>10000) break;


    }

    EvTree->Write();
 	MyFile->Close();
 	cout << "number MC " << MC_aboveQ2 << endl;
 	cout << "number analyzed " << track << endl;
 	cout << "Track 3 events: " << alltrack3 << endl;
 	cout << "Mass cut NOT accepted: " << notaccepted << endl;
	
}





#include "jpsi.h"
#include "BranchesInfo.h"
#include "caloheader.h"
#include "clusterizer.cxx"

void jpsi(
	//TString inFile            	= "../macros/18x272-IP6-Bill.root",
	TString inFile            	= "./../rootfiles/all.root",
	TString inFileGeometry      = "./../rootfiles/geometry.root",
	bool do_reclus              = true,
    unsigned short primaryTrackSource = 0,
    bool HOFrame  				= true

){


	// load tree
    TChain *const tt_event = new TChain("T");
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
    tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("geometry_tree");
    if(!tt_geometry){ cout << "geometry tree not found... returning!"<< endl; return;}

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    std::cout << "Number of events in tree: " << nEntriesTree << std::endl;

    SetBranchAddressesTree(tt_event);
    SetBranchAddressesGeometryTree(tt_geometry);
    SetGeometryIndices();


    TFile *MyFile = new TFile("./../rootfiles/jpsi_final.root","RECREATE");
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

	double p_eta_HO, p_theta_HO, p_p_HO; 
	double p_pT_HO, p_E_HO, p_phi_HO; 

	double _Q2,_xbj, _phiJpsi, _phiJpsiMC;
	double _t, _Q2MC, _xbjMC, _tMC, _y, _yMC;
	double _xv, _xvMC;
	double _missmass, _missmassp;

	int _BECAL_nclusters;

	float* _BECAL_E   = new float[_BECAL_nclusters];
	float* _BECAL_Eta = new float[_BECAL_nclusters];
	float* _BECAL_Phi = new float[_BECAL_nclusters];
	float* _BECAL_x   = new float[_BECAL_nclusters];
	float* _BECAL_y   = new float[_BECAL_nclusters];
	float* _BECAL_z   = new float[_BECAL_nclusters];
	int* _BECAL_ID  = new int[_BECAL_nclusters];
	int* _BECAL_NtrueID =  new int[_BECAL_nclusters];

	int _FEMC_nclusters;

	float* _FEMC_E   = new float[_FEMC_nclusters];
	float* _FEMC_Eta = new float[_FEMC_nclusters];
	float* _FEMC_Phi = new float[_FEMC_nclusters];
	float* _FEMC_x   = new float[_FEMC_nclusters];
	float* _FEMC_y   = new float[_FEMC_nclusters];
	float* _FEMC_z   = new float[_FEMC_nclusters];
	int* _FEMC_ID  = new int[_FEMC_nclusters];
	int* _FEMC_NtrueID =  new int[_FEMC_nclusters];

	int _EEMC_nclusters;

	float* _EEMC_E   = new float[_EEMC_nclusters];
	float* _EEMC_Eta = new float[_EEMC_nclusters];
	float* _EEMC_Phi = new float[_EEMC_nclusters];
	float* _EEMC_x   = new float[_EEMC_nclusters];
	float* _EEMC_y   = new float[_EEMC_nclusters];
	float* _EEMC_z   = new float[_EEMC_nclusters];
	int* _EEMC_ID  = new int[_EEMC_nclusters];
	int* _EEMC_NtrueID =  new int[_EEMC_nclusters];


	int nHits; 
	int hits_layerID[20];
	int hits_trueID[20];

	EvTree->Branch("Jpsi_track3_M1",&Jpsi_track3_M1,"Jpsi_track3_M1/D");
	EvTree->Branch("Jpsi_track3_M2",&Jpsi_track3_M2,"Jpsi_track3_M2/D"); 	
	
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
	

	EvTree->Branch("RP1",&_RP1,"RP1/I");
    EvTree->Branch("RP2",&_RP2,"RP2/I");
    EvTree->Branch("RPhits",&_RPhits,"RPhits/I");
    EvTree->Branch("RPx",_RPx,"RPx[RPhits]/F");
    EvTree->Branch("RPy",_RPy,"RPy[RPhits]/F");
    EvTree->Branch("RPz",_RPz,"RPz[RPhits]/F");
   	EvTree->Branch("RPind",_RPind,"RPind[RPhits]/I");
    EvTree->Branch("RPtrPx",_RPtrPx,"RPtrPx[RPhits]/F");
   	EvTree->Branch("RPtrPy",_RPtrPy,"RPtrPy[RPhits]/F");
    EvTree->Branch("RPtrPz",_RPtrPz,"RPtrPz[RPhits]/F");
    EvTree->Branch("RPid",_RPid,"RPid[RPhits]/I");
    EvTree->Branch("RPpx",_RPpx,"RPpx[RPhits]/F");
    EvTree->Branch("RPpy",_RPpy,"RPpy[RPhits]/F");
    EvTree->Branch("RPpz",_RPpz,"RPpz[RPhits]/F");
    EvTree->Branch("RPpid",_RPpid,"RPpid[RPhits]/I");

    EvTree->Branch("B0hits", &_B0hits, "B0hits/I");
    EvTree->Branch("B0x", _B0x, "B0x[B0hits]/F");
    EvTree->Branch("B0y", _B0y, "B0y[B0hits]/F");
    EvTree->Branch("B0z", _B0z, "B0z[B0hits]/F");
    EvTree->Branch("B0ind", _B0ind, "B0ind[B0hits]/I");
    EvTree->Branch("B0px", _B0px, "B0px[B0hits]/F");
    EvTree->Branch("B0py", _B0py, "B0py[B0hits]/F");
    EvTree->Branch("B0pz", _B0pz, "B0pz[B0hits]/F");
    EvTree->Branch("B0trPx", _B0trPx, "B0trPx[B0hits]/F");
    EvTree->Branch("B0trPy", _B0trPy, "B0trPy[B0hits]/F");
    EvTree->Branch("B0trPz", _B0trPz, "B0trPz[B0hits]/F");
    EvTree->Branch("B0id", _B0id, "B0id[B0hits]/I");

	EvTree->Branch("Q2",&_Q2,"Q2/D");
	EvTree->Branch("Q2MC",&_Q2MC,"Q2MC/D");
	EvTree->Branch("t",&_t,"t/D");
	EvTree->Branch("tMC",&_tMC,"tMC/D");
	EvTree->Branch("xbj",&_xbj,"xbj/D");
	EvTree->Branch("xbjMC",&_xbjMC,"xbjMC/D");
	EvTree->Branch("phiJpsi",&_phiJpsi,"phiJpsi/D");
	EvTree->Branch("phiJpsiMC",&_phiJpsiMC,"phiJpsiMC/D"); 
    EvTree->Branch("y",&_y,"y/D");
    EvTree->Branch("yMC",&_yMC,"yMC/D");
    EvTree->Branch("xv",&_xv,"xv/D");
	EvTree->Branch("xvMC",&_xvMC,"xvMC/D");


    EvTree->Branch("nTracks",&_nTracks,"nTracks/I");
    EvTree->Branch("BECAL_nclusters", &_BECAL_nclusters, "BECAL_nclusters/I");
    EvTree->Branch("BECAL_E ", _BECAL_E, "BECAL_E[BECAL_nclusters]/F");
    EvTree->Branch("BECAL_Eta", _BECAL_Eta, "BECAL_Eta[BECAL_nclusters]/F");
    EvTree->Branch("BECAL_Phi", _BECAL_Phi, "BECAL_Phi[BECAL_nclusters]/F");
    EvTree->Branch("BECAL_x", _BECAL_x, "BECAL_x[BECAL_nclusters]/F");
    EvTree->Branch("BECAL_y", _BECAL_y, "BECAL_y[BECAL_nclusters]/F");
    EvTree->Branch("BECAL_z", _BECAL_z, "BECAL_z[BECAL_nclusters]/F");
    EvTree->Branch("BECAL_ID", _BECAL_ID, "BECAL_ID[BECAL_nclusters]/I");
    EvTree->Branch("BECAL_NtrueID", _BECAL_NtrueID, "BECAL_NtrueID[BECAL_nclusters]/I");

	EvTree->Branch("EEMC_nclusters", &_EEMC_nclusters, "EEMC_nclusters/I");
    EvTree->Branch("EEMC_E ", _EEMC_E, "EEMC_E[EEMC_nclusters]/F");
    EvTree->Branch("EEMC_Eta", _EEMC_Eta, "EEMC_Eta[EEMC_nclusters]/F");
    EvTree->Branch("EEMC_Phi", _EEMC_Phi, "EEMC_Phi[EEMC_nclusters]/F");
    EvTree->Branch("EEMC_x", _EEMC_x, "EEMC_x[EEMC_nclusters]/F");
    EvTree->Branch("EEMC_y", _EEMC_y, "EEMC_y[EEMC_nclusters]/F");
    EvTree->Branch("EEMC_z", _EEMC_z, "EEMC_z[EEMC_nclusters]/F");
    EvTree->Branch("EEMC_ID", _EEMC_ID, "EEMC_ID[EEMC_nclusters]/I");
    EvTree->Branch("EEMC_NtrueID", _EEMC_NtrueID, "EEMC_NtrueID[EEMC_nclusters]/I");

	EvTree->Branch("FEMC_nclusters", &_FEMC_nclusters, "FEMC_nclusters/I");
    EvTree->Branch("FEMC_E ", _FEMC_E, "FEMC_E[FEMC_nclusters]/F");
    EvTree->Branch("FEMC_Eta", _FEMC_Eta, "FEMC_Eta[FEMC_nclusters]/F");
    EvTree->Branch("FEMC_Phi", _FEMC_Phi, "FEMC_Phi[FEMC_nclusters]/F");
    EvTree->Branch("FEMC_x", _FEMC_x, "FEMC_x[FEMC_nclusters]/F");
    EvTree->Branch("FEMC_y", _FEMC_y, "FEMC_y[FEMC_nclusters]/F");
    EvTree->Branch("FEMC_z", _FEMC_z, "FEMC_z[FEMC_nclusters]/F");
    EvTree->Branch("FEMC_ID", _FEMC_ID, "FEMC_ID[FEMC_nclusters]/I");
    EvTree->Branch("FEMC_NtrueID", _FEMC_NtrueID, "FEMC_NtrueID[FEMC_nclusters]/I");
    
    EvTree->Branch("missmass", &_missmass, "missmass/D");
    EvTree->Branch("missmassp", &_missmassp, "missmassp/D");


	if(HepmcEnabled){
    		AddBranchesHepmc(EvTree);
	}



   	_nEventsTree=0;

    int track = 0;
    //int MC_aboveQ2 = 0;

    TLorentzVector ProtonBeam(0,0,274.998,275);
    TLorentzVector ElectronBeam(0,0,-18,18);
    TLorentzRotation rotlor = TLorentzRotation().RotateY(12.5e-3).Boost(sin(12.5e-3),0,0);

    int notaccepted =0 ;
    int alltrack3 = 0 ;

    // main event loop
//    for (Long64_t i=0; i<1000; i++) {  
    for (Long64_t i=0; i<nEntriesTree;i++) {  

    	tt_event->GetEntry(i);
    	
    	if (_hepmcp_Q2 < 1.) continue;
    	
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
    		}

    		else continue;


    	}
    	
    

		TLorentzVector hypep;
		TLorentzVector hyp1;
		TLorentzVector hyp2;

		int eptrack_N[3];
		int othercharged[3];
		int noe = 0;
		int nep = 0;
		//if(_nTracks != 3) continue;

		alltrack3++;
		
		if(_nTracks ==3){
			
	
			//======= Read from the tracking paths ==============//
			for(int l=0; l<_nTracks; l++){
				//if(_track_source[l]!=0) continue;
		
				if( _track_charge[l] == short(1)) {
					eptrack_N[nep] = l;
					nep++; 

				}
				else if( _track_charge[l] == short(-1))  {
					othercharged[noe] = l;
					noe++;
				}
			}
		

		if(nep>1) continue;
		if(noe!=2) continue;
		
		
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



	float seed_E = 0.5;
        float aggregation_E = 0.1;

        EEMC_cluster.clear();
        FHCAL_cluster.clear();
        BECAL_cluster.clear();
        HCALIN_cluster.clear();
        HCALOUT_cluster.clear();
        FEMC_cluster.clear();
        EHCAL_cluster.clear();


        float seed_E_EEMC = 0.1;
        float aggregation_E_EEMC = 0.05;
		float seed_E_FHCAL = 0.5;
        float aggregation_E_FHCAL = 0.1;
        float seed_E_FEMC         = 0.1;
        float aggregation_E_FEMC  = 0.005;
        float seed_E_HCALIN = 0.2;
        float aggregation_E_HCALIN = 0.05;
        float seed_E_HCALOUT = 0.5;
        float aggregation_E_HCALOUT = 0.1;
        float seed_E_BECAL = 0.1;
        float aggregation_E_BECAL = 0.01;
        float seed_E_EHCAL = 0.01;
        float aggregation_E_EHCAL = 0.005;

	runclusterizer(kMA, kEEMC,    seed_E_EEMC, 	aggregation_E_EEMC, 0);
	runclusterizer(kMA, kFEMC,    seed_E_FEMC, aggregation_E_FEMC, 0);
	runclusterizer(kMA, kBECAL,   seed_E_BECAL, aggregation_E_BECAL, 0);
		
	
	runclusterizer(kMA, kFHCAL,   seed_E_FHCAL, aggregation_E_FHCAL, 0);
	runclusterizer(kMA, kHCALIN,  seed_E_HCALIN, aggregation_E_HCALIN, 0);
	runclusterizer(kMA, kHCALOUT, seed_E_HCALOUT, aggregation_E_HCALOUT, 0);
	runclusterizer(kMA, kEHCAL,   seed_E_EHCAL, aggregation_E_EHCAL, 0);

	_BECAL_nclusters = 0;
	map<int, clustersStrct>::iterator itr;
	for (itr = BECAL_cluster.begin(); itr != BECAL_cluster.end(); ++itr) {
		if(itr->second.cluster_E > emin){
		_BECAL_E[_BECAL_nclusters] =  itr->second.cluster_E;
		_BECAL_Eta[_BECAL_nclusters] = itr->second.cluster_Eta;
		_BECAL_Phi[_BECAL_nclusters] = itr->second.cluster_Phi;
		_BECAL_x[_BECAL_nclusters] = itr->second.cluster_X;
		_BECAL_y[_BECAL_nclusters] = itr->second.cluster_Y;
		_BECAL_z[_BECAL_nclusters] = itr->second.cluster_Z;
		_BECAL_ID[_BECAL_nclusters] = itr->second.cluster_trueID;
		_BECAL_NtrueID[_BECAL_nclusters] = itr->second.cluster_NtrueID;
		_BECAL_nclusters++;
		}		
	}

	_EEMC_nclusters = 0;
	map<int, clustersStrct>::iterator itr2;
	for (itr2 = EEMC_cluster.begin(); itr2 != EEMC_cluster.end(); ++itr2) {
		if(itr2->second.cluster_E > emin){
		_EEMC_E[_EEMC_nclusters] =  itr2->second.cluster_E;
		_EEMC_Eta[_EEMC_nclusters] = itr2->second.cluster_Eta;
		_EEMC_Phi[_EEMC_nclusters] = itr2->second.cluster_Phi;
		_EEMC_x[_EEMC_nclusters] = itr2->second.cluster_X;
		_EEMC_y[_EEMC_nclusters] = itr2->second.cluster_Y;
		_EEMC_z[_EEMC_nclusters] = itr2->second.cluster_Z;
		_EEMC_ID[_EEMC_nclusters] = itr2->second.cluster_trueID;
		_EEMC_NtrueID[_EEMC_nclusters] = itr2->second.cluster_NtrueID;
		_EEMC_nclusters++;	
		}	
	}

	_FEMC_nclusters = 0;
	map<int, clustersStrct>::iterator itr3;
	for (itr3 = FEMC_cluster.begin(); itr3 != FEMC_cluster.end(); ++itr3) {
		if(itr3->second.cluster_E > emin){
		_FEMC_E[_FEMC_nclusters] =  itr3->second.cluster_E;
		_FEMC_Eta[_FEMC_nclusters] = itr3->second.cluster_Eta;
		_FEMC_Phi[_FEMC_nclusters] = itr3->second.cluster_Phi;
		_FEMC_x[_FEMC_nclusters] = itr3->second.cluster_X;
		_FEMC_y[_FEMC_nclusters] = itr3->second.cluster_Y;
		_FEMC_z[_FEMC_nclusters] = itr3->second.cluster_Z;
		_FEMC_ID[_FEMC_nclusters] = itr3->second.cluster_trueID;
		_FEMC_NtrueID[_FEMC_nclusters] = itr3->second.cluster_NtrueID;
		_FEMC_nclusters++;		
	}
	}
		double res_track = 0.02;
		if(_EEMC_nclusters>1){
			for(int n=0; n<_EEMC_nclusters; n++)

				if( fabs( (_EEMC_Eta[n] - Electron.Eta() )/Electron.Eta()) < res_track &&  fabs((_EEMC_Phi[n] - Electron.Phi())/Electron.Phi() )  < res_track && _EEMC_Eta[n] < -2.  && _EEMC_Eta[n] > -3.5 ){
				
        	      	double pp 	 = sqrt(_EEMC_E[n]*_EEMC_E[n] - me*me);
					double theta = 2*atan(exp(-_EEMC_Eta[n]));

					double px = pp*sin(theta)*cos(_EEMC_Phi[n]);
					double py = pp*sin(theta)*sin(_EEMC_Phi[n]);
					double pz = pp*cos(theta);					

         			Electron.SetPxPyPzE(px, py, pz, _EEMC_E[n]);

         			//cout << Electron.M() << endl;
         			break;

				}
		}




		//======================= Proton reconstructed ===============


		TLorentzVector proton;
		TVector3 prot3D(_RPtrPx[0],_RPtrPy[0],_RPtrPz[0]);
		//TVector3 prot3D(_RPpx[0],_RPpy[0],_RPpz[0]);
		double  energyproton = sqrt(prot3D.Mag2() + protonmass*protonmass );
		proton.SetPxPyPzE(_RPtrPx[0], _RPtrPy[0], _RPtrPz[0], energyproton);
		//proton.SetPxPyPzE(_RPpx[0], _RPpy[0], _RPpz[0], energyproton);




		if(HOFrame){
			JPsi = rotlor*JPsi;
			em   = rotlor*em;
			ep   = rotlor*ep;
			Electron = rotlor*Electron;
			proton = rotlor*proton;
		}

		
		Jpsi_M = JPsi.M();
		ep_eta = hypep.Eta(); em_eta = em.Eta(); Jpsi_eta = JPsi.Eta(); e_eta = Electron.Eta();
		ep_theta = hypep.Theta(); em_theta = em.Theta(); Jpsi_theta = JPsi.Theta(); e_theta = Electron.Theta();
		ep_phi = hypep.Phi(); em_phi = em.Phi(); Jpsi_phi = JPsi.Phi(); e_phi = Electron.Phi();
		ep_p = hypep.Vect().Mag(); em_p = em.Vect().Mag(); Jpsi_p = JPsi.Vect().Mag(); e_p = Electron.Vect().Mag();
		ep_pT = hypep.Perp(); em_pT = em.Perp(); Jpsi_pT = JPsi.Perp(); e_pT = Electron.Perp();
		ep_E = hypep.E(); em_E = em.E(); Jpsi_E = JPsi.E(); e_E = Electron.E();
		p_eta = proton.Eta(); p_theta = proton.Theta(); p_phi = proton.Phi(); p_p = proton.Vect().Mag();
		p_pT = proton.Perp(); p_E = proton.E(); 


		//==== Kinematics from tracking ==========

	   	


		_Q2MC = -(ElectronBeam - eeHepmc)*(ElectronBeam - eeHepmc);
		_Q2   = -(ElectronBeam - Electron)*(ElectronBeam - Electron);


		_xbjMC = _Q2MC/(2*ProtonBeam*(ElectronBeam - eeHepmc));
		_xbj   = _Q2/(2*ProtonBeam*(ElectronBeam - Electron));

		_tMC = (pHepmc - ProtonBeam)*(pHepmc - ProtonBeam);
		_t = (proton - ProtonBeam)*(proton - ProtonBeam);

		//======phi 
		TVector3 q 		  = (ElectronBeam - Electron).Vect();
		TVector3 v1      =  q.Cross(Electron.Vect());
		TVector3 v2      =  q.Cross(JPsi.Vect());
		_phiJpsi          =  v1.Angle(v2);

		if(q.Dot(v1.Cross(v2))<0) _phiJpsi=2.*TMath::Pi()-_phiJpsi;

		TVector3 qHepmc		  = (ElectronBeam - eeHepmc).Vect();
		TVector3 v1Hepmc      =  q.Cross(eeHepmc.Vect());
		TVector3 v2Hepmc      =  q.Cross(JpsiHepmc.Vect());
		_phiJpsiMC          =  v1.Angle(v2);

		if(q.Dot(v1.Cross(v2))<0) _phiJpsiMC=2.*TMath::Pi()-_phiJpsiMC;


		//====== Rapidity 

		_y    = log((proton.E() + proton.Pz())/(proton.E() - proton.Pz()))/2;
		_yMC  = log((pHepmc.E() + pHepmc.Pz())/(pHepmc.E() - pHepmc.Pz()))/2;


		//====== Rapidity 
		_xv     = (_Q2MC + JPsi.M())/(2*ProtonBeam*(ElectronBeam - eeHepmc));
		_xvMC   = (_Q2  + JpsiHepmc.M())/(2*ProtonBeam*(ElectronBeam - Electron));

		_missmass  = (ElectronBeam +  ProtonBeam - JPsi - proton - Electron).M();

	    _missmassp  = (ElectronBeam +  ProtonBeam - JPsi  - Electron).M();
	}

	 


		track++;
	

	  	EvTree->Fill();
    	//if (i>10000) break;


    }

    EvTree->Write();
 	MyFile->Close();
 	//cout << "number MC " << MC_aboveQ2 << endl;
 	cout << "number analyzed " << track << endl;
 	cout << "Track 3 events: " << alltrack3 << endl;
 	cout << "Mass cut NOT accepted: " << notaccepted << endl;
	
}


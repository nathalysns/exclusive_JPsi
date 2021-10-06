#include <TROOT.h>
#include <TString.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TVector3.h>


float idep = 7;
float idem = 8; 
const double mp =  0.93827208816;    
const double me =  0.000511;
const double low_limjpsiM = 2.;
const double high_limjpsiM = 5.;


const int _maxNTowers = 50 * 50;
const int _maxNTowersCentral = 2000;
const int _maxNTowersDR = 3000 * 3000;
const int _maxNTowersCalo = 5000000;
const int _maxNclusters = 100;
const int _maxNclustersCentral = 2000;
const int _maxNTracks = 200;
const int _maxNProjections = 2000;
const int _maxNMCPart = 100000;
const int _maxNHepmcp = 1000;
const int _maxNHits = 1000;

float _nEventsTree;

enum calotype {
  kFHCAL          = 0,
  kFEMC           = 1,
  kDRCALO         = 2,
  kEEMC           = 3,
  kCEMC           = 4,
  kEHCAL          = 5,
  kHCALIN         = 6,
  kHCALOUT        = 7,
  kLFHCAL         = 8,
  kEEMCG          = 9,
  kBECAL          = 10
};

bool caloEnabled[20]      = {0};
bool tracksEnabled        = 0;
bool vertexEnabled        = 0;
bool xSectionEnabled      = 0;
bool HepmcEnabled         = 0;
bool HitsEnabled          = 0;

// Event level info
float _cross_section;
float _event_weight;
int _n_generator_accepted;

// track hits
//int _nHitsLayers;
//int* _hits_layerID              = new int[_maxNHits];
//float* _hits_x             = new float[_maxNHits];
//float* _hits_y             = new float[_maxNHits];
//float* _hits_z             = new float[_maxNHits];
//float* _hits_t             = new float[_maxNHits];

// tracks
int _nTracks;
float* _track_ID                 = new float[_maxNTracks];
float* _track_trueID             = new float[_maxNTracks];
float* _track_px                 = new float[_maxNTracks];
float* _track_py                 = new float[_maxNTracks];
float* _track_pz                 = new float[_maxNTracks];
short* _track_charge             = new short[_maxNTracks];
unsigned short* _track_source             = new unsigned short[_maxNTracks];

// vertex
float _vertex_x;
float _vertex_y;
float _vertex_z;

int _nProjections;
float* _track_ProjTrackID                 = new float[_maxNProjections];
int* _track_ProjLayer                 = new int[_maxNProjections];
float* _track_Proj_x             = new float[_maxNProjections];
float* _track_Proj_y             = new float[_maxNProjections];
float* _track_Proj_z             = new float[_maxNProjections];
float* _track_Proj_t             = new float[_maxNProjections];
float* _track_Proj_true_x             = new float[_maxNProjections];
float* _track_Proj_true_y             = new float[_maxNProjections];
float* _track_Proj_true_z             = new float[_maxNProjections];
float* _track_Proj_true_t             = new float[_maxNProjections];
float* _tracks_source                 = new float[_maxNProjections];

// MC particles
int _nMCPart;
int* _mcpart_ID                  = new int[_maxNMCPart];
int* _mcpart_ID_parent           = new int[_maxNMCPart];
int* _mcpart_PDG                 = new int[_maxNMCPart];
float* _mcpart_E                 = new float[_maxNMCPart];
float* _mcpart_px                = new float[_maxNMCPart];
float* _mcpart_py                = new float[_maxNMCPart];
float* _mcpart_pz                = new float[_maxNMCPart];
float* _mcpart_Eta                = new float[_maxNMCPart];
float* _mcpart_Phi                = new float[_maxNMCPart];
int* _mcpart_BCID                  = new int[_maxNMCPart];

// towers
int _nTowers_CEMC;
float* _tower_CEMC_E            = new float[_maxNTowersCentral];
int* _tower_CEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_CEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_CEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EEMC;
float* _tower_EEMC_E            = new float[_maxNTowersCentral];
int* _tower_EEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_EEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_EEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EEMCG;
float* _tower_EEMCG_E            = new float[_maxNTowersCentral];
int* _tower_EEMCG_iEta         = new int[_maxNTowersCentral];
int* _tower_EEMCG_iPhi         = new int[_maxNTowersCentral];
int* _tower_EEMCG_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EHCAL;
float* _tower_EHCAL_E            = new float[_maxNTowersCentral];
int* _tower_EHCAL_iEta         = new int[_maxNTowersCentral];
int* _tower_EHCAL_iPhi         = new int[_maxNTowersCentral];
int* _tower_EHCAL_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALIN;
float* _tower_HCALIN_E            = new float[_maxNTowersCentral];
int* _tower_HCALIN_iEta         = new int[_maxNTowersCentral];
int* _tower_HCALIN_iPhi         = new int[_maxNTowersCentral];
int* _tower_HCALIN_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALOUT;
float* _tower_HCALOUT_E            = new float[_maxNTowersCentral];
int* _tower_HCALOUT_iEta         = new int[_maxNTowersCentral];
int* _tower_HCALOUT_iPhi         = new int[_maxNTowersCentral];
int* _tower_HCALOUT_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_DRCALO;
float* _tower_DRCALO_E            = new float[_maxNTowersDR];
int* _tower_DRCALO_iEta         = new int[_maxNTowersDR];
int* _tower_DRCALO_iPhi         = new int[_maxNTowersDR];
int* _tower_DRCALO_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_LFHCAL;
float* _tower_LFHCAL_E            = new float[_maxNTowersDR];
int* _tower_LFHCAL_iEta         = new int[_maxNTowersDR];
int* _tower_LFHCAL_iPhi         = new int[_maxNTowersDR];
int* _tower_LFHCAL_iL           = new int[_maxNTowersDR];
int* _tower_LFHCAL_trueID       = new int[_maxNTowersDR];


// towers
int _nTowers_BECAL;
float* _tower_BECAL_E            = new float[_maxNTowersDR];
int* _tower_BECAL_iEta         = new int[_maxNTowersDR];
int* _tower_BECAL_iPhi         = new int[_maxNTowersDR];
int* _tower_BECAL_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_FHCAL;
float* _tower_FHCAL_E            = new float[_maxNTowers];
int* _tower_FHCAL_iEta         = new int[_maxNTowers];
int* _tower_FHCAL_iPhi         = new int[_maxNTowers];
int* _tower_FHCAL_trueID       = new int[_maxNTowers];

// towers
int _nTowers_FEMC;
float* _tower_FEMC_E             = new float[_maxNTowers];
int* _tower_FEMC_iEta          = new int[_maxNTowers];
int* _tower_FEMC_iPhi          = new int[_maxNTowers];
int* _tower_FEMC_trueID        = new int[_maxNTowers];

// clusters
int _nclusters_FHCAL;
float* _clusters_FHCAL_E            = new float[_maxNclusters];
float* _clusters_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_FHCAL_Phi         = new float[_maxNclusters];
int* _clusters_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_FHCAL_trueID       = new int[_maxNclusters];

// clusters
int _nclusters_FEMC;
float* _clusters_FEMC_E             = new float[_maxNclusters];
float* _clusters_FEMC_Eta          = new float[_maxNclusters];
float* _clusters_FEMC_Phi          = new float[_maxNclusters];
int* _clusters_FEMC_NTower          = new int[_maxNclusters];
int* _clusters_FEMC_trueID        = new int[_maxNclusters];

// clusters
int _nclusters_BECAL;
float* _clusters_BECAL_E            = new float[_maxNclusters];
float* _clusters_BECAL_Eta         = new float[_maxNclusters];
float* _clusters_BECAL_Phi         = new float[_maxNclusters];
int* _clusters_BECAL_NTower         = new int[_maxNclusters];
int* _clusters_BECAL_trueID       = new int[_maxNclusters];

int _nHepmcp; 
int _hepmcp_procid; 
float _hepmcp_x1; 
float _hepmcp_x2;
float _hepmcp_Q2;
int* _hepmcp_status            = new int[_maxNHepmcp];
int* _hepmcp_PDG               = new int[_maxNHepmcp];
float* _hepmcp_E               = new float[_maxNHepmcp]; 
float* _hepmcp_px              = new float[_maxNHepmcp]; 
float* _hepmcp_py              = new float[_maxNHepmcp];  
float* _hepmcp_pz              = new float[_maxNHepmcp];
int* _hepmcp_BCID              = new int[_maxNHepmcp]; 
int* _hepmcp_m1                = new int[_maxNHepmcp];
int* _hepmcp_m2                = new int[_maxNHepmcp]; 

int _nHits;
int *_hits_layerID     = new int[_maxNHits];
int *_hits_trueID      = new int[_maxNHits];
float* _hits_x          = new float[_maxNHits];
float* _hits_y          = new float[_maxNHits];
float* _hits_z          = new float[_maxNHits];
float* _hits_x2         = new float[_maxNHits];
float* _hits_y2         = new float[_maxNHits];
float* _hits_z2         = new float[_maxNHits];
float* _hits_t          = new float[_maxNHits];
float* _hits_edep       = new float[_maxNHits];

void SetBranchAddressesTree(TTree* inputTree){

    if (inputTree->GetBranchStatus("cross_section") ){
      xSectionEnabled = 1;
      inputTree->SetBranchAddress("cross_section", &_cross_section);
      inputTree->SetBranchAddress("event_weight", &_event_weight);
      inputTree->SetBranchAddress("n_generator_accepted", &_n_generator_accepted);
    }
    /*inputTree->SetBranchAddress("nHits",                        &_nHitsLayers);
    inputTree->SetBranchAddress("hits_layerID",                 _hits_layerID);
    inputTree->SetBranchAddress("hits_x",               _hits_x);
    inputTree->SetBranchAddress("hits_y",               _hits_y);
    inputTree->SetBranchAddress("hits_z",               _hits_z);
    inputTree->SetBranchAddress("hits_t",               _hits_t);
	*/
    if (inputTree->GetBranchStatus("nTracks") ){
      tracksEnabled = 1;
      inputTree->SetBranchAddress("nTracks",              &_nTracks);
      inputTree->SetBranchAddress("track_ID",            _track_ID);
      inputTree->SetBranchAddress("track_px",            _track_px);
      inputTree->SetBranchAddress("track_py",            _track_py);
      inputTree->SetBranchAddress("track_pz",            _track_pz);
      inputTree->SetBranchAddress("track_trueID",        _track_trueID);
     // inputTree->SetBranchAddress("track_source",        _track_source);
       inputTree->SetBranchAddress("tracks_charge",       _track_charge);
/*
      inputTree->SetBranchAddress("nProjections",         &_nProjections);
      inputTree->SetBranchAddress("track_ProjTrackID",    _track_ProjTrackID);
      inputTree->SetBranchAddress("track_ProjLayer",      _track_ProjLayer);

      inputTree->SetBranchAddress("track_TLP_x",           _track_Proj_x);
      inputTree->SetBranchAddress("track_TLP_y",           _track_Proj_y);
      inputTree->SetBranchAddress("track_TLP_z",           _track_Proj_z);
      //inputTree->SetBranchAddress("track_TLP_t",           _track_Proj_t);
      inputTree->SetBranchAddress("track_TLP_true_x",      _track_Proj_true_x);
      inputTree->SetBranchAddress("track_TLP_true_y",      _track_Proj_true_y);
      inputTree->SetBranchAddress("track_TLP_true_z",      _track_Proj_true_z);
      //inputTree->SetBranchAddress("track_TLP_true_t",      _track_Proj_true_t);
  */    
    }

    // towers EEMC
    if( inputTree->GetBranchStatus("tower_EEMC_N") ){
      caloEnabled[kEEMC] = 1;
      inputTree->SetBranchAddress("tower_EEMC_N",                &_nTowers_EEMC);
      inputTree->SetBranchAddress("tower_EEMC_E",                _tower_EEMC_E);
      inputTree->SetBranchAddress("tower_EEMC_iEta",             _tower_EEMC_iEta);
      inputTree->SetBranchAddress("tower_EEMC_iPhi",             _tower_EEMC_iPhi);
      inputTree->SetBranchAddress("tower_EEMC_trueID",           _tower_EEMC_trueID);
    } 
    // towers EEMCG
    if( inputTree->GetBranchStatus("tower_EEMCG_N") ){
      caloEnabled[kEEMCG] = 1;
      inputTree->SetBranchAddress("tower_EEMCG_N",                &_nTowers_EEMCG);
      inputTree->SetBranchAddress("tower_EEMCG_E",                _tower_EEMCG_E);
      inputTree->SetBranchAddress("tower_EEMCG_iEta",             _tower_EEMCG_iEta);
      inputTree->SetBranchAddress("tower_EEMCG_iPhi",             _tower_EEMCG_iPhi);
      inputTree->SetBranchAddress("tower_EEMCG_trueID",           _tower_EEMCG_trueID);
    } 

    // towers EHCAL
    if( inputTree->GetBranchStatus("tower_EHCAL_N") ){
      caloEnabled[kEHCAL] = 1;
      inputTree->SetBranchAddress("tower_EHCAL_N",                &_nTowers_EHCAL);
      inputTree->SetBranchAddress("tower_EHCAL_E",                _tower_EHCAL_E);
      inputTree->SetBranchAddress("tower_EHCAL_iEta",             _tower_EHCAL_iEta);
      inputTree->SetBranchAddress("tower_EHCAL_iPhi",             _tower_EHCAL_iPhi);
      inputTree->SetBranchAddress("tower_EHCAL_trueID",           _tower_EHCAL_trueID);
    }

    // towers HCALIN
    if( inputTree->GetBranchStatus("tower_HCALIN_N") ){
      caloEnabled[kHCALIN] = 1;
      inputTree->SetBranchAddress("tower_HCALIN_N",                &_nTowers_HCALIN);
      inputTree->SetBranchAddress("tower_HCALIN_E",                _tower_HCALIN_E);
      inputTree->SetBranchAddress("tower_HCALIN_iEta",             _tower_HCALIN_iEta);
      inputTree->SetBranchAddress("tower_HCALIN_iPhi",             _tower_HCALIN_iPhi);
      inputTree->SetBranchAddress("tower_HCALIN_trueID",           _tower_HCALIN_trueID);
    } 

    // towers HCALOUT
    if( inputTree->GetBranchStatus("tower_HCALOUT_N") ){
      caloEnabled[kHCALOUT] = 1;
      inputTree->SetBranchAddress("tower_HCALOUT_N",                &_nTowers_HCALOUT);
      inputTree->SetBranchAddress("tower_HCALOUT_E",                _tower_HCALOUT_E);
      inputTree->SetBranchAddress("tower_HCALOUT_iEta",             _tower_HCALOUT_iEta);
      inputTree->SetBranchAddress("tower_HCALOUT_iPhi",             _tower_HCALOUT_iPhi);
      inputTree->SetBranchAddress("tower_HCALOUT_trueID",           _tower_HCALOUT_trueID);
    } 

	// towers CEMC
    if( inputTree->GetBranchStatus("tower_CEMC_N") ){
      caloEnabled[kCEMC] = 1;
      inputTree->SetBranchAddress("tower_CEMC_N",                &_nTowers_CEMC);
      inputTree->SetBranchAddress("tower_CEMC_E",                _tower_CEMC_E);
      inputTree->SetBranchAddress("tower_CEMC_iEta",             _tower_CEMC_iEta);
      inputTree->SetBranchAddress("tower_CEMC_iPhi",             _tower_CEMC_iPhi);
      inputTree->SetBranchAddress("tower_CEMC_trueID",           _tower_CEMC_trueID);
    } 

    // towers BECAL
    if( inputTree->GetBranchStatus("tower_BECAL_N") ){
      caloEnabled[kBECAL] = 1;
      inputTree->SetBranchAddress("tower_BECAL_N",                &_nTowers_BECAL);
      inputTree->SetBranchAddress("tower_BECAL_E",                _tower_BECAL_E);
      inputTree->SetBranchAddress("tower_BECAL_iEta",             _tower_BECAL_iEta);
      inputTree->SetBranchAddress("tower_BECAL_iPhi",             _tower_BECAL_iPhi);
      inputTree->SetBranchAddress("tower_BECAL_trueID",           _tower_BECAL_trueID);
    } 

    // towers LFHCAL
    if( inputTree->GetBranchStatus("tower_LFHCAL_N") ){
      caloEnabled[kLFHCAL] = 1;
      inputTree->SetBranchAddress("tower_LFHCAL_N",                &_nTowers_LFHCAL);
      inputTree->SetBranchAddress("tower_LFHCAL_E",                _tower_LFHCAL_E);
      inputTree->SetBranchAddress("tower_LFHCAL_iEta",             _tower_LFHCAL_iEta);
      inputTree->SetBranchAddress("tower_LFHCAL_iPhi",             _tower_LFHCAL_iPhi);
      inputTree->SetBranchAddress("tower_LFHCAL_iL",               _tower_LFHCAL_iL);
      inputTree->SetBranchAddress("tower_LFHCAL_trueID",           _tower_LFHCAL_trueID);
    } 
    // towers DRCALO
    if( inputTree->GetBranchStatus("tower_DRCALO_N") ){
      caloEnabled[kDRCALO] = 1;
      inputTree->SetBranchAddress("tower_DRCALO_N",                &_nTowers_DRCALO);
      inputTree->SetBranchAddress("tower_DRCALO_E",                _tower_DRCALO_E);
      inputTree->SetBranchAddress("tower_DRCALO_iEta",             _tower_DRCALO_iEta);
      inputTree->SetBranchAddress("tower_DRCALO_iPhi",             _tower_DRCALO_iPhi);
      inputTree->SetBranchAddress("tower_DRCALO_trueID",           _tower_DRCALO_trueID);
    } 

    // towers FHCAL
    if( inputTree->GetBranchStatus("tower_FHCAL_N") ){
      caloEnabled[kFHCAL] = 1;
      inputTree->SetBranchAddress("tower_FHCAL_N",                &_nTowers_FHCAL);
      inputTree->SetBranchAddress("tower_FHCAL_E",                _tower_FHCAL_E);
      inputTree->SetBranchAddress("tower_FHCAL_iEta",             _tower_FHCAL_iEta);
      inputTree->SetBranchAddress("tower_FHCAL_iPhi",             _tower_FHCAL_iPhi);
      inputTree->SetBranchAddress("tower_FHCAL_trueID",           _tower_FHCAL_trueID);
    } 

    // towers FEMC
    if( inputTree->GetBranchStatus("tower_FEMC_N") ){
      caloEnabled[kFEMC] = 1;
      inputTree->SetBranchAddress("tower_FEMC_N",                 &_nTowers_FEMC);
      inputTree->SetBranchAddress("tower_FEMC_E",                 _tower_FEMC_E);
      inputTree->SetBranchAddress("tower_FEMC_iEta",              _tower_FEMC_iEta);
      inputTree->SetBranchAddress("tower_FEMC_iPhi",              _tower_FEMC_iPhi);
      inputTree->SetBranchAddress("tower_FEMC_trueID",            _tower_FEMC_trueID);
    }

    // clusters HCAL
    if (caloEnabled[kFHCAL]){
      inputTree->SetBranchAddress("cluster_FHCAL_N",                &_nclusters_FHCAL);
      inputTree->SetBranchAddress("cluster_FHCAL_E",                _clusters_FHCAL_E);
      inputTree->SetBranchAddress("cluster_FHCAL_Eta",             _clusters_FHCAL_Eta);
      inputTree->SetBranchAddress("cluster_FHCAL_Phi",             _clusters_FHCAL_Phi);
      inputTree->SetBranchAddress("cluster_FHCAL_NTower",             _clusters_FHCAL_NTower);
      inputTree->SetBranchAddress("cluster_FHCAL_trueID",           _clusters_FHCAL_trueID);
    }
    // clusters EMC
    if (caloEnabled[kFEMC]){
      inputTree->SetBranchAddress("cluster_FEMC_N",                 &_nclusters_FEMC);
      inputTree->SetBranchAddress("cluster_FEMC_E",                 _clusters_FEMC_E);
      inputTree->SetBranchAddress("cluster_FEMC_Eta",              _clusters_FEMC_Eta);
      inputTree->SetBranchAddress("cluster_FEMC_Phi",              _clusters_FEMC_Phi);
      inputTree->SetBranchAddress("cluster_FEMC_NTower",              _clusters_FEMC_NTower);
      inputTree->SetBranchAddress("cluster_FEMC_trueID",            _clusters_FEMC_trueID);
    }
    // vertex

    if (inputTree->GetBranchStatus("vertex_x") ){
      vertexEnabled = 1;
      inputTree->SetBranchAddress("vertex_x",                     &_vertex_x);
      inputTree->SetBranchAddress("vertex_y",                     &_vertex_y);
      inputTree->SetBranchAddress("vertex_z",                     &_vertex_z);
    }
    // MC particles
    inputTree->SetBranchAddress("nMCPart",       &_nMCPart);
    inputTree->SetBranchAddress("mcpart_ID",     _mcpart_ID);
    inputTree->SetBranchAddress("mcpart_ID_parent",     _mcpart_ID_parent);
    inputTree->SetBranchAddress("mcpart_PDG",    _mcpart_PDG);
    inputTree->SetBranchAddress("mcpart_E",      _mcpart_E);
    inputTree->SetBranchAddress("mcpart_px",     _mcpart_px);
    inputTree->SetBranchAddress("mcpart_py",     _mcpart_py);
    inputTree->SetBranchAddress("mcpart_pz",     _mcpart_pz);
    inputTree->SetBranchAddress("mcpart_BCID",     _mcpart_BCID);

    if (inputTree->GetBranchStatus("maxNHepmcp") ){
      HepmcEnabled = 1;
      inputTree->SetBranchAddress("maxNHepmcp",      &_nHepmcp);
      //inputTree->SetBranchAddress("hepmcp_procid",      &_hepmcp_procid);
      inputTree->SetBranchAddress("hepmcp_x1",          &_hepmcp_x1);
      inputTree->SetBranchAddress("hepmcp_x2",          &_hepmcp_x2);
      inputTree->SetBranchAddress("hepmcp_Q2",          &_hepmcp_Q2);
      inputTree->SetBranchAddress("hepmcp_status",      _hepmcp_status);
      inputTree->SetBranchAddress("hepmcp_PDG",         _hepmcp_PDG);
      inputTree->SetBranchAddress("hepmcp_E",           _hepmcp_E);
      inputTree->SetBranchAddress("hepmcp_px",          _hepmcp_px);
      inputTree->SetBranchAddress("hepmcp_py",          _hepmcp_py);
      inputTree->SetBranchAddress("hepmcp_pz",          _hepmcp_pz);
      inputTree->SetBranchAddress("hepmcp_BCID",        _hepmcp_BCID);
      inputTree->SetBranchAddress("hepmcp_m1",          _hepmcp_m1);
      inputTree->SetBranchAddress("hepmcp_m2",          _hepmcp_m2);
    }

    if (inputTree->GetBranchStatus("nHits") ){
      HitsEnabled = 1;
      inputTree->SetBranchAddress("nHits",      &_nHits);
      inputTree->SetBranchAddress("hits_layerID",   _hits_layerID);
      inputTree->SetBranchAddress("hits_trueID",    _hits_trueID);
      inputTree->SetBranchAddress("hits_x",         _hits_x);
      inputTree->SetBranchAddress("hits_y",         _hits_y);
      inputTree->SetBranchAddress("hits_z",         _hits_z);
      inputTree->SetBranchAddress("hits_x2",        _hits_x2);
      inputTree->SetBranchAddress("hits_y2",        _hits_y2);
      inputTree->SetBranchAddress("hits_z2",        _hits_z2);
      inputTree->SetBranchAddress("hits_t",         _hits_t);
      inputTree->SetBranchAddress("hits_edep",      _hits_edep);
    }
}

TRandom3  _fRandom;                                  // random for effi generation

TTree * tt_geometry;
int _calogeom_ID;
int _calogeom_towers_N;
int*  _calogeom_towers_iEta = new int[_maxNTowersCalo];
int*  _calogeom_towers_iPhi = new int[_maxNTowersCalo];
int*  _calogeom_towers_iL   = new int[_maxNTowersCalo];
float*  _calogeom_towers_Eta = new float[_maxNTowersCalo];
float*  _calogeom_towers_Phi = new float[_maxNTowersCalo];
float*  _calogeom_towers_x = new float[_maxNTowersCalo];
float*  _calogeom_towers_y = new float[_maxNTowersCalo];
float*  _calogeom_towers_z = new float[_maxNTowersCalo];

void SetBranchAddressesGeometryTree(TTree* inputTreeGeo){
    inputTreeGeo->SetBranchAddress("calo",              &_calogeom_ID);
    inputTreeGeo->SetBranchAddress("calo_towers_N",     &_calogeom_towers_N);
    inputTreeGeo->SetBranchAddress("calo_towers_iEta",  _calogeom_towers_iEta);
    inputTreeGeo->SetBranchAddress("calo_towers_iPhi",  _calogeom_towers_iPhi);
    inputTreeGeo->SetBranchAddress("calo_towers_iL",    _calogeom_towers_iL);
    inputTreeGeo->SetBranchAddress("calo_towers_Eta",   _calogeom_towers_Eta);
    inputTreeGeo->SetBranchAddress("calo_towers_Phi",   _calogeom_towers_Phi);
    inputTreeGeo->SetBranchAddress("calo_towers_x",     _calogeom_towers_x);
    inputTreeGeo->SetBranchAddress("calo_towers_y",     _calogeom_towers_y);
    inputTreeGeo->SetBranchAddress("calo_towers_z",     _calogeom_towers_z);
}

int GetCorrectMCArrayEntry(float objectTrueID){
  for(Int_t imc=0; imc<_nMCPart; imc++){
    if(objectTrueID==_mcpart_ID[imc]){
      return imc;
    }
  }
  return -1;
}

TLorentzVector rotate_reco(TLorentzVector init, double angle ){
  TLorentzVector rotat;
  TLorentzRotation l = TLorentzRotation().RotateY(angle/2).Boost(sin(angle/2),0,0);
  rotat = l*init;
  return rotat;
}



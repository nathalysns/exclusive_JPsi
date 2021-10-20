//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in excl_ana.h.
//
// excl_ana(const std::string &name = "excl_ana")
// everything is keyed to excl_ana, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// excl_ana::~excl_ana()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int excl_ana::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int excl_ana::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int excl_ana::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int excl_ana::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int excl_ana::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int excl_ana::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int excl_ana::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void excl_ana::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "excl_ana.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <g4eval/CaloEvalStack.h>
//#include <calobase/RawCluster.h>
#include <g4eval/CaloRawClusterEval.h>
//#include <calobase/RawClusterContainer.h>
#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>

#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


/// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <g4eval/SvtxEvalStack.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>


/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include "g4eval/CaloRawTowerEval.h"

using namespace std;


//____________________________________________________________________________..
//excl_ana::excl_ana(const std::string &name):
// SubsysReco(name)
//{
//  std::cout << "excl_ana::excl_ana(const std::string &name) Calling ctor" << std::endl;
//}


excl_ana::excl_ana(const std::string &name):
 SubsysReco(name)
 , _caloevalstackFHCAL(nullptr)
  , _caloevalstackBECAL(nullptr)
  , _caloevalstackHCALIN(nullptr)
  , _caloevalstackHCALOUT(nullptr)
  , _caloevalstackEHCAL(nullptr)
  , _caloevalstackDRCALO(nullptr)
  , _caloevalstackFOCAL(nullptr)
  , _caloevalstackLFHCAL(nullptr)
  , _caloevalstackFEMC(nullptr)
  , _caloevalstackCEMC(nullptr)
  , _caloevalstackEEMC(nullptr)
  , _caloevalstackEEMCG(nullptr)
 , _strict(false)
{

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

  set_do_PROJECTIONS(false);
  initializeVariables();
  initializeTrees();

  _reco_e_threshold[kFHCAL]   = 0.05;
  _reco_e_threshold[kFEMC]    = 0.005;
  _reco_e_threshold[kDRCALO]  = 0.0;
  _reco_e_threshold[kFOCAL]  = 0.0;
  _reco_e_threshold[kEEMC]    = 0.005;
  _reco_e_threshold[kCEMC]    = 0.01;
  _reco_e_threshold[kEHCAL]   = 0.05;
  _reco_e_threshold[kHCALIN]  = 0.01;
  _reco_e_threshold[kHCALOUT] = 0.05;
  _reco_e_threshold[kLFHCAL]  = 0.001;
  _reco_e_threshold[kEEMCG]   = 0.005;
  _reco_e_threshold[kBECAL]   = 0.001;

}


void excl_ana::SetOutputFile(TString outfilename) {

        mOutputFileName = outfilename;

}


//____________________________________________________________________________..
excl_ana::~excl_ana()
{
  delete tree;

if(_svtxEvalStack) delete _svtxEvalStack;
  gsl_rng_free(m_RandomGenerator);

  std::cout << "excl_ana::~excl_ana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int excl_ana::Init(PHCompositeNode *topNode)
{

  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(mOutputFileName, "RECREATE");

  std::cout << "excl_ana::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::InitRun(PHCompositeNode *topNode)
{
  std::cout << "excl_ana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::process_event(PHCompositeNode *topNode)
{
  
//===================== Number of events read
   event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;




//=========================
/*
_nMCPart = 0;
PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");

if (truthinfocontainer){
   PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetParticleRange();

    for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr){
     PHG4Particle* g4particle = truth_itr->second;
     if (!g4particle) continue;
     if (g4particle->get_e() < 2) continue;
     if (  g4particle->get_pid() <-12 && g4particle->get_pid() > 2500) continue;
    	_mcpart_ID[_nMCPart] = g4particle->get_track_id();
     	_mcpart_ID_parent[_nMCPart] = g4particle->get_parent_id();
     	_mcpart_PDG[_nMCPart] = g4particle->get_pid();
     	_mcpart_E[_nMCPart] = g4particle->get_e();
     	_mcpart_px[_nMCPart] = g4particle->get_px();
     	_mcpart_py[_nMCPart] = g4particle->get_py();
     	_mcpart_pz[_nMCPart] = g4particle->get_pz();

     	_mcpart_BCID[_nMCPart] = g4particle->get_barcode();

     	_nMCPart++;
    }
}else{
     	cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
}*/

// ========= lopping over HEPMC Particles =============================================//

PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

 _nHepmcp = 0 ;
if (hepmceventmap){
  if (Verbosity() > 0) cout << "saving HepMC output" << endl;

  for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin(); eventIter != hepmceventmap->end(); ++eventIter){
           PHHepMCGenEvent* hepmcevent = eventIter->second;
    if (hepmcevent){
      HepMC::GenEvent* truthevent = hepmcevent->getEvent();
      if (!truthevent){
        cout << PHWHERE
        << "no evt pointer under phhepmvgeneventmap found "
 	<< endl;
 	return 0;
      }

     HepMC::PdfInfo* pdfinfo = truthevent->pdf_info();

      _hepmcp_x1 = pdfinfo->x1();
      _hepmcp_x2 = pdfinfo->x2();
      _hepmcp_Q2 = pdfinfo->scalePDF();
 
      for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
 	   iter != truthevent->particles_end(); ++iter ){

 	  _hepmcp_E[_nHepmcp] = (*iter)->momentum().e();
      	  _hepmcp_PDG[_nHepmcp] = (*iter)->pdg_id();
  	  _hepmcp_px[_nHepmcp] = (*iter)->momentum().px();
          _hepmcp_py[_nHepmcp] = (*iter)->momentum().py();
          _hepmcp_pz[_nHepmcp] = (*iter)->momentum().pz();
          _hepmcp_status[_nHepmcp] = (*iter)->status();
          _hepmcp_BCID[_nHepmcp] = (*iter)->barcode();
          _hepmcp_m2[_nHepmcp] = 0;
          _hepmcp_m1[_nHepmcp] = 0;

          if ((*iter)->production_vertex()){
               
	     for (HepMC::GenVertex::particle_iterator mother = (*iter)->production_vertex()->particles_begin(HepMC::parents);
                  mother != (*iter)->production_vertex()->particles_end(HepMC::parents);++mother){
                 _hepmcp_m2[_nHepmcp] = (*mother)->barcode();
                
		 if (_hepmcp_m1[_nHepmcp] == 0)
                   _hepmcp_m1[_nHepmcp] = (*mother)->barcode();
               	 }
           }//end production vertex 
           if (Verbosity() > 2) cout << "nHepmcp " << _nHepmcp << "\tPDG " << _hepmcp_PDG[_nHepmcp] << "\tEnergy " << _hepmcp_E[_nHepmcp] << "\tbarcode " << _hepmcp_BCID[_nHepmcp] << "\tMother1 " << _hepmcp_m1[_nHepmcp]<< "\tMother2 " << _hepmcp_m2[_nHepmcp] << endl;
           _nHepmcp++;
        }//end for Hparticle_const_iterator
     }// if hepmcevent
  }//end eventIter loo 
  _maxNHepmcp = _nHepmcp;
}//end if loop 
else {
  if (Verbosity() > 0) cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
  return 0;
}


//========== trackinf ====================================//
_nTracks = 0;
_nProjections = 0;  
// Loop over track maps, identifiy each source.
// Although this configuration is fixed here, it doesn't require multiple sources.// It will only store them if they're available.

std::vector<std::pair<std::string, TrackSource_t>> trackMapInfo = {{"TrackMap", TrackSource_t::all}, {"TrackMapInner", TrackSource_t::inner}}; 

bool foundAtLeastOneTrackSource = false;
for (const auto& trackMapInfo : trackMapInfo){
   SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, trackMapInfo.first);
   if (trackmap){
      foundAtLeastOneTrackSource = true;
      int nTracksInASource = 0;
      if (Verbosity() > 0) cout << "saving tracks for track map: " << trackMapInfo.first << endl;
      for (SvtxTrackMap::ConstIter track_itr = trackmap->begin(); track_itr != trackmap->end(); track_itr++){
         SvtxTrack_FastSim* track = dynamic_cast<SvtxTrack_FastSim*>(track_itr->second);
         if (track){
             _track_ID[_nTracks] = track->get_id();
	    
             _track_px[_nTracks] = track->get_px();
             _track_py[_nTracks] = track->get_py();
             _track_pz[_nTracks] = track->get_pz();
             _track_charge[_nTracks] = track->get_charge();
             // Ideally, would be dca3d_xy and dca3d_z, but these don't seem to be calculated properly in the        
              // current (June 2021) simulations (they return NaN). So we take dca (seems to be ~ the 3d distance)
              // and dca_2d (seems to be ~ the distance in the transverse plane).
              // The names of the branches are based on the method names.
             _track_dca[_nTracks] = static_cast<float>(track->get_dca());
	     _track_dca_2d[_nTracks] = static_cast<float>(track->get_dca2d());
             _track_trueID[_nTracks] = track->get_truth_track_id();
             _track_source[_nTracks] = static_cast<unsigned short>(trackMapInfo.second);
	     if (_do_PROJECTIONS){
		// find projections
		for (SvtxTrack::ConstStateIter trkstates = track->begin_states(); trkstates != track->end_states(); ++trkstates){
                   if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " processing " << trkstates->second->get_name() << endl;
	           string trackStateName = trkstates->second->get_name();
		   if (Verbosity() > 1)   cout << __PRETTY_FUNCTION__ << " found " << trkstates->second->get_name() << endl;
		   int trackStateIndex = GetProjectionIndex(trackStateName);
		   if (trackStateIndex > -1){
	 	       // save true projection info to given branch
	 	       _track_TLP_true_x[_nProjections] = trkstates->second->get_pos(0);
		       _track_TLP_true_y[_nProjections] = trkstates->second->get_pos(1);
                       _track_TLP_true_z[_nProjections] = trkstates->second->get_pos(2);
		       _track_TLP_true_t[_nProjections] = trkstates->first;
		       _track_ProjLayer[_nProjections] = trackStateIndex;
                       _track_ProjTrackID[_nProjections] = _nTracks;
    
                       string nodename = "G4HIT_" + trkstates->second->get_name();
                       PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
		       if (hits){
			   if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
			   PHG4HitContainer::ConstRange hit_range = hits->getHits();
			   for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
			      if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;
			      if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0){
			         if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
		                 // save reco projection info to given branch
		                 _track_TLP_x[_nProjections] = hit_iter->second->get_x(0);
	                         _track_TLP_y[_nProjections] = hit_iter->second->get_y(0);
			         _track_TLP_z[_nProjections] = hit_iter->second->get_z(0);
                                _track_TLP_t[_nProjections] = hit_iter->second->get_t(0);
			      }// if get_trkid
			    }// for hit iter
		        } //if hits
			else {
			  if (Verbosity() > 1) if (Verbosity() > 1);
			  continue;
			}
			_nProjections++;
                    } //trackstand indesx		
	         }//trackstates
	      }// do projections
	      _nTracks++;
	     nTracksInASource++;
           } else {

	        if (Verbosity() > 0) {
		   cout << "PHG4TrackFastSimEval::fill_track_tree - ignore track that is not a SvtxTrack_FastSim:";
		  track_itr->second->identify();
		}
		continue;
	   }
	   if (Verbosity() > 0) cout << "saved\t" << nTracksInASource << "\ttracks from track map " << trackMapInfo.first << ". Total saved tracks: " << _nTracks << endl;

	}// end track iter
}         else {

	if (Verbosity() > 0) cout << PHWHERE << "SvtxTrackMap node with name '" << trackMapInfo.first << "' not found on node tree" << endl;
	}
 //  }//end if trackmap
    if (foundAtLeastOneTrackSource == false) {
        cout << PHWHERE << "Requested tracks, but found no sources on node tree. Returning" << endl;
        return 0;
    }
}//end for tracMapinfo
 
 _maxNProjections = _nProjections;
 _maxNTracks = _nTracks;

//===================================================== B0

 process_B0(topNode);
 process_RomanPots(topNode,1);
 process_ZDC(topNode);

  //----------------------
  //    TOWERS FHCAL
  //----------------------


  CaloEvalStack* _caloevalstackFHCAL;
  _caloevalstackFHCAL = new CaloEvalStack(topNode, "FHCAL");
  _caloevalstackFHCAL->set_strict(_strict);
  _caloevalstackFHCAL->set_verbosity(Verbosity() + 1);

  CaloRawTowerEval* towerevalFHCAL = _caloevalstackFHCAL->get_rawtower_eval();
  _nTowers_FHCAL = 0;
 
  for (Int_t itow = 0; itow < _maxNTowers; itow++)
  {
    _tower_FHCAL_E[itow] = 0;
    _tower_FHCAL_iEta[itow] = 0;
    _tower_FHCAL_iPhi[itow] = 0;
    _tower_FHCAL_trueID[itow] = 0;
  }

  string towernodeFHCAL = "TOWER_CALIB_FHCAL";
  RawTowerContainer* towersFHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeFHCAL.c_str());
  if (towersFHCAL)
  {
    if (Verbosity() > 0)
    {
      cout << "saving HCAL towers" << endl;
    }
    string towergeomnodeFHCAL = "TOWERGEOM_FHCAL";
    RawTowerGeomContainer* towergeomFHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeFHCAL.c_str());
    if (towergeomFHCAL)
    {
        RawTowerContainer::ConstRange begin_end = towersFHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            if (tower->get_energy() < _reco_e_threshold[kFHCAL]) continue;

            _tower_FHCAL_iEta[_nTowers_FHCAL] = tower->get_bineta();
            _tower_FHCAL_iPhi[_nTowers_FHCAL] = tower->get_binphi();
            _tower_FHCAL_E[_nTowers_FHCAL] = tower->get_energy();
	
	    PHG4Particle* primary = towerevalFHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_FHCAL_trueID[_nTowers_FHCAL] = primary->get_track_id();
	     }
            else
            {
              _tower_FHCAL_trueID[_nTowers_FHCAL] = -10;
            }
            _nTowers_FHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeFHCAL << endl;
        }
        }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_FHCAL << "\tFHCAL towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeFHCAL << endl;
      }
  } 

  //----------------------
  //  TOWERS BECAL
  //----------------------

    if (!_caloevalstackBECAL)
    {
      _caloevalstackBECAL = new CaloEvalStack(topNode, "BECAL");
      _caloevalstackBECAL->set_strict(_strict);
      _caloevalstackBECAL->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackBECAL->next_event(topNode);
    }
   //========= initial
    _nTowers_BECAL = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_BECAL_E[itow] = 0;
      _tower_BECAL_iEta[itow] = 0;
      _tower_BECAL_iPhi[itow] = 0;
      _tower_BECAL_trueID[itow] = 0;
    }


    CaloRawTowerEval* towerevalBECAL = _caloevalstackBECAL->get_rawtower_eval();
    _nTowers_BECAL = 0;
    string towernodeBECAL = "TOWER_CALIB_BECAL";
    RawTowerContainer* towersBECAL = findNode::getClass<RawTowerContainer>(topNode, towernodeBECAL.c_str());
    
    if (Verbosity() > 0)
    {
      RawTowerContainer* towersBECAL1 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_RAW_BECAL");
      RawTowerContainer* towersBECAL2 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_SIM_BECAL");
      cout << "BECAL sim: " << towersBECAL2->size() << endl;
      cout << "BECAL raw: " << towersBECAL1->size() << endl;
      cout << "BECAL calib: " << towersBECAL->size() << endl;
    }
    if (towersBECAL)
    {
      if (Verbosity() > 0)
      {
        cout << "saving BECAL towers" << endl;
      }
      string towergeomnodeBECAL = "TOWERGEOM_BECAL";
      RawTowerGeomContainer* towergeomBECAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeBECAL.c_str());
      if (towergeomBECAL)
      {
       RawTowerContainer::ConstRange begin_end = towersBECAL->getTowers();
       RawTowerContainer::ConstIterator rtiter;
       for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
       {
          RawTower* tower = rtiter->second;
          if (tower)
          {
	     // min energy cut
	      if (tower->get_energy() < _reco_e_threshold[kBECAL]) continue;
            _tower_BECAL_iEta[_nTowers_BECAL] = tower->get_bineta();
            _tower_BECAL_iPhi[_nTowers_BECAL] = tower->get_binphi();
            _tower_BECAL_E[_nTowers_BECAL] = tower->get_energy();

            PHG4Particle* primary = towerevalBECAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_BECAL_trueID[_nTowers_BECAL] = primary->get_track_id();
            }
            else
            {
              _tower_BECAL_trueID[_nTowers_BECAL] = -10;
            }
            _nTowers_BECAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeBECAL << endl;
        }
     }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_BECAL << "\tBECAL towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeBECAL << endl;
      }
   }

  //----------------------
  //    TOWERS HCALIN
  //----------------------
   _nTowers_HCALIN = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_HCALIN_E[itow] = 0;
      _tower_HCALIN_iEta[itow] = 0;
      _tower_HCALIN_iPhi[itow] = 0;
      _tower_HCALIN_trueID[itow] = 0;
    }
    if (!_caloevalstackHCALIN)
    {
      _caloevalstackHCALIN = new CaloEvalStack(topNode, "HCALIN");
      _caloevalstackHCALIN->set_strict(_strict);
      _caloevalstackHCALIN->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackHCALIN->next_event(topNode);
    } 
  
   CaloRawTowerEval* towerevalHCALIN = _caloevalstackHCALIN->get_rawtower_eval();
   string towernodeHCALIN = "TOWER_CALIB_HCALIN";
   RawTowerContainer* towersHCALIN = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALIN.c_str());
   if (towersHCALIN)
    {
      if (Verbosity() > 0)
      {
        cout << "saving HCAL towers" << endl;
      }
      string towergeomnodeHCALIN = "TOWERGEOM_HCALIN";
      RawTowerGeomContainer* towergeomHCALIN = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALIN.c_str());
      if (towergeomHCALIN)
      {
	RawTowerContainer::ConstRange begin_end = towersHCALIN->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
	    // min energy cut
	    if (tower->get_energy() < _reco_e_threshold[kHCALIN]) continue;
            _tower_HCALIN_iEta[_nTowers_HCALIN] = tower->get_bineta();
            _tower_HCALIN_iPhi[_nTowers_HCALIN] = tower->get_binphi();
            _tower_HCALIN_E[_nTowers_HCALIN] = tower->get_energy();
	
	    PHG4Particle* primary = towerevalHCALIN->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_HCALIN_trueID[_nTowers_HCALIN] = primary->get_track_id();
	    }
            else
            {
              _tower_HCALIN_trueID[_nTowers_HCALIN] = -10;
            }
            _nTowers_HCALIN++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALIN << endl;
        }
       }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_HCALIN << "\tHCALIN towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeHCALIN << endl;
      }
   }	
   //----------------------
   //    TOWERS HCALOUT
   //----------------------
   
    _nTowers_HCALOUT = 0;
    string towernodeHCALOUT = "TOWER_CALIB_HCALOUT";
    RawTowerContainer* towersHCALOUT = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALOUT.c_str());

    if (!_caloevalstackHCALOUT)
    {
      _caloevalstackHCALOUT = new CaloEvalStack(topNode, "HCALOUT");
      _caloevalstackHCALOUT->set_strict(_strict);
      _caloevalstackHCALOUT->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackHCALOUT->next_event(topNode);
    }

    if (towersHCALOUT)
    {
      if (Verbosity() > 0)
      {
        cout << "saving HCAL towers" << endl;
      }
      string towergeomnodeHCALOUT = "TOWERGEOM_HCALOUT";
      RawTowerGeomContainer* towergeomHCALOUT = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALOUT.c_str());
      CaloRawTowerEval* towerevalHCALOUT = _caloevalstackHCALOUT->get_rawtower_eval();
      if (towergeomHCALOUT)
      {
        RawTowerContainer::ConstRange begin_end = towersHCALOUT->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
 	   // min energy cut
 	   if (tower->get_energy() < _reco_e_threshold[kHCALOUT]) continue;
            _tower_HCALOUT_iEta[_nTowers_HCALOUT] = tower->get_bineta();
            _tower_HCALOUT_iPhi[_nTowers_HCALOUT] = tower->get_binphi();
            _tower_HCALOUT_E[_nTowers_HCALOUT] = tower->get_energy();

	    PHG4Particle* primary = towerevalHCALOUT->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = primary->get_track_id();
	    }
            else
            {
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = -10;
            }

            _nTowers_HCALOUT++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALOUT << endl;
        }
       }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_HCALOUT << "\tHCALOUT towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeHCALOUT << endl;
      }
   }

  //----------------------
  //    TOWERS EHCAL
  //----------------------

    _nTowers_EHCAL = 0;
    if (!_caloevalstackEHCAL)
    {
      _caloevalstackEHCAL = new CaloEvalStack(topNode, "EHCAL");
      _caloevalstackEHCAL->set_strict(_strict);
      _caloevalstackEHCAL->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackEHCAL->next_event(topNode);
    }

    string towernodeEHCAL = "TOWER_CALIB_EHCAL";
    RawTowerContainer* towersEHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeEHCAL.c_str());
    if (towersEHCAL)
    {
      if (Verbosity() > 0)
      {
        cout << "saving HCAL towers" << endl;
      }
      string towergeomnodeEHCAL = "TOWERGEOM_EHCAL";
      RawTowerGeomContainer* towergeomEHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEHCAL.c_str());
      CaloRawTowerEval* towerevalEHCAL = _caloevalstackEHCAL->get_rawtower_eval();
      if (towergeomEHCAL)
      {
	RawTowerContainer::ConstRange begin_end = towersEHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
	    // min energy cut
	    if (tower->get_energy() < _reco_e_threshold[kEHCAL]) continue;
            _tower_EHCAL_iEta[_nTowers_EHCAL] = tower->get_bineta();
            _tower_EHCAL_iPhi[_nTowers_EHCAL] = tower->get_binphi();
            _tower_EHCAL_E[_nTowers_EHCAL] = tower->get_energy();

	    PHG4Particle* primary = towerevalEHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_EHCAL_trueID[_nTowers_EHCAL] = primary->get_track_id();
	    }
            else
            {
              _tower_EHCAL_trueID[_nTowers_EHCAL] = -10;
            }

            _nTowers_EHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEHCAL << endl;
        }
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_EHCAL << "\tEHCAL towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeEHCAL << endl;
      }

    }

    _nTowers_LFHCAL = 0;
    if (!_caloevalstackLFHCAL)
    {
      _caloevalstackLFHCAL = new CaloEvalStack(topNode, "LFHCAL");
      _caloevalstackLFHCAL->set_strict(_strict);
      _caloevalstackLFHCAL->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackLFHCAL->next_event(topNode);
    }
    string towernodeLFHCAL = "TOWER_CALIB_LFHCAL";
    RawTowerContainer* towersLFHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeLFHCAL.c_str());
    if (Verbosity() > 1) 
      std::cout << "reading towers: "<< towersLFHCAL->size() << std::endl;
    if (towersLFHCAL)
    {
      if (Verbosity() > 0)
      {
        cout << "saving LFHCAL towers" << endl;
      }
      string towergeomnodeLFHCAL = "TOWERGEOM_LFHCAL";
      RawTowerGeomContainer* towergeomLFHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeLFHCAL.c_str());
      CaloRawTowerEval* towerevalLFHCAL = _caloevalstackLFHCAL->get_rawtower_eval();
      if (towergeomLFHCAL)
      {
         if (Verbosity() > 0)
        {
          cout << "found LFHCAL geom" << endl;
        }

        RawTowerContainer::ConstRange begin_end = towersLFHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
	    if (tower->get_energy() < _reco_e_threshold[kLFHCAL]) continue; 
            if (Verbosity() > 1) cout << "\n event eval: \t" << tower->get_energy()<< "\t ieta: " << tower->get_bineta()<< "\t iphi: " << tower->get_binphi() << "\t iZ: " << tower->get_binl()<< endl;
            _tower_LFHCAL_iEta[_nTowers_LFHCAL] = tower->get_bineta();
            _tower_LFHCAL_iPhi[_nTowers_LFHCAL] = tower->get_binphi();
            _tower_LFHCAL_iL[_nTowers_LFHCAL] = tower->get_binl();
            _tower_LFHCAL_E[_nTowers_LFHCAL] = tower->get_energy();

	    PHG4Particle* primary = towerevalLFHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_LFHCAL_trueID[_nTowers_LFHCAL] = primary->get_track_id();
	    }
            else
            {
              _tower_LFHCAL_trueID[_nTowers_LFHCAL] = -10;
            }

            _nTowers_LFHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeLFHCAL << endl;
        }
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_LFHCAL << "\tLFHCAL towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeLFHCAL << endl;
      }
    }
   //----------------------
   //    TOWERS FEMC
   //---------------------

    _nTowers_FEMC = 0;
    if (!_caloevalstackFEMC)
    {
      _caloevalstackFEMC = new CaloEvalStack(topNode, "FEMC");
      _caloevalstackFEMC->set_strict(_strict);
      _caloevalstackFEMC->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackFEMC->next_event(topNode);
    }
    string towernodeFEMC = "TOWER_CALIB_FEMC";
    RawTowerContainer* towersFEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeFEMC.c_str());
    if (towersFEMC)
    {
      if (Verbosity() > 0)
      {
        cout << "saving EMC towers" << endl;
      }
      string towergeomnodeFEMC = "TOWERGEOM_FEMC";
      CaloRawTowerEval* towerevalFEMC = _caloevalstackFEMC->get_rawtower_eval();
      RawTowerContainer::ConstRange begin_end = towersFEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
	    // min energy cut
	    if (tower->get_energy() < _reco_e_threshold[kFEMC]) continue;

            _tower_FEMC_iEta[_nTowers_FEMC] = tower->get_bineta();
            _tower_FEMC_iPhi[_nTowers_FEMC] = tower->get_binphi();
            _tower_FEMC_E[_nTowers_FEMC] = tower->get_energy();
  	     PHG4Particle* primary = towerevalFEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_FEMC_trueID[_nTowers_FEMC] = primary->get_track_id();
	     }
            else
            {
              _tower_FEMC_trueID[_nTowers_FEMC] = -10;
            }
	    	    
	     _nTowers_FEMC++;
          }
        }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeFEMC << endl;
      }
    }
   
  //----------------------
  //    TOWERS EEMC
  //----------------------

    _nTowers_EEMC = 0;
    if (!_caloevalstackEEMC)
    {
      _caloevalstackEEMC = new CaloEvalStack(topNode, "EEMC");
      _caloevalstackEEMC->set_strict(_strict);
      _caloevalstackEEMC->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackEEMC->next_event(topNode);
    }
    string towernodeEEMC = "TOWER_CALIB_EEMC";
    RawTowerContainer* towersEEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeEEMC.c_str());
    if (towersEEMC)
    {
      if (Verbosity() > 0)
      {
        cout << "saving EMC towers" << endl;
      }
      string towergeomnodeEEMC = "TOWERGEOM_EEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEEMC.c_str());
      CaloRawTowerEval* towerevalEEMC = _caloevalstackEEMC->get_rawtower_eval();
      if (towergeom)
      {
      RawTowerContainer::ConstRange begin_end = towersEEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            if (tower->get_energy() < _reco_e_threshold[kEEMC]) continue;

            _tower_EEMC_iEta[_nTowers_EEMC] = tower->get_bineta();
            _tower_EEMC_iPhi[_nTowers_EEMC] = tower->get_binphi();
            _tower_EEMC_E[_nTowers_EEMC] = tower->get_energy();
	    PHG4Particle* primary = towerevalEEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_EEMC_trueID[_nTowers_EEMC] = primary->get_track_id();
	    }
            else
            {
              _tower_EEMC_trueID[_nTowers_EEMC] = -10;
            }
            _nTowers_EEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0) 
        {
	   cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEEMC << endl;
        }
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_EEMC << "\tEEMC towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeEEMC << endl;
      }
    }
   
  //----------------------
  //    TOWERS EEMC
  //----------------------
  
  _nTowers_EEMCG = 0;
  if (!_caloevalstackEEMCG)
  {
    _caloevalstackEEMCG = new CaloEvalStack(topNode, "EEMC_glass");
    _caloevalstackEEMCG->set_strict(_strict);
    _caloevalstackEEMCG->set_verbosity(Verbosity() + 1);
  }
  else
  {
    _caloevalstackEEMCG->next_event(topNode);
  }
  string towernodeEEMCG = "TOWER_CALIB_EEMC_glass";
  RawTowerContainer* towersEEMCG = findNode::getClass<RawTowerContainer>(topNode, towernodeEEMCG.c_str());
  if (towersEEMCG)
  {
    if (Verbosity() > 0)
    {
      cout << "saving EMC towers" << endl;
    }
    string towergeomnodeEEMCG = "TOWERGEOM_EEMC_glass";
    RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEEMCG.c_str());
    CaloRawTowerEval* towerevalEEMCG = _caloevalstackEEMCG->get_rawtower_eval();
    if (towergeom)
    {
	RawTowerContainer::ConstRange begin_end = towersEEMCG->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
	   // min energy cut
	  if (tower->get_energy() < _reco_e_threshold[kEEMCG]) continue;

           _tower_EEMCG_iEta[_nTowers_EEMCG] = tower->get_bineta();
           _tower_EEMCG_iPhi[_nTowers_EEMCG] = tower->get_binphi();
           _tower_EEMCG_E[_nTowers_EEMCG] = tower->get_energy();
	   PHG4Particle* primary = towerevalEEMCG->max_truth_primary_particle_by_energy(tower);
           if (primary)
           {
             _tower_EEMCG_trueID[_nTowers_EEMCG] = primary->get_track_id();
	   }
            else
            {
              _tower_EEMCG_trueID[_nTowers_EEMCG] = -10;
            }
           _nTowers_EEMCG++;
         }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEEMCG << endl;
        }
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_EEMCG << "\tEEMCG towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeEEMCG << endl;
      }
   }


  tree->Fill();


//cout<<" Event filled, next "<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "excl_ana::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::EndRun(const int runnumber)
{
  std::cout << "excl_ana::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::End(PHCompositeNode *topNode)
{
  std::cout << "excl_ana::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  
  outfile->cd();
  tree->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  if (_caloevalstackFHCAL) delete _caloevalstackFHCAL;
  if (_caloevalstackBECAL) delete _caloevalstackBECAL;
  if (_caloevalstackHCALIN) delete _caloevalstackHCALIN;
  if (_caloevalstackHCALOUT) delete _caloevalstackHCALOUT;
  if (_caloevalstackEHCAL) delete _caloevalstackEHCAL;
  if (_caloevalstackDRCALO) delete _caloevalstackDRCALO;
  if (_caloevalstackFOCAL) delete _caloevalstackFOCAL;
  if (_caloevalstackLFHCAL) delete _caloevalstackLFHCAL;
  if (_caloevalstackFEMC) delete _caloevalstackFEMC;
  if (_caloevalstackCEMC) delete _caloevalstackCEMC;
  if (_caloevalstackEEMC) delete _caloevalstackEEMC;
  if (_caloevalstackEEMCG) delete _caloevalstackEEMCG;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int excl_ana::Reset(PHCompositeNode *topNode)
{
 std::cout << "excl_ana::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void excl_ana::Print(const std::string &what) const
{
  std::cout << "excl_ana::Print(const std::string &what) const Printing info for " << what << std::endl;
}


//***********************************************************
void excl_ana::initializeTrees()
{
  tree = new TTree("T", "A tree for DVCS");

  tree->Branch("nHepmcp", &_nHepmcp, "nHepmcp/I");
  tree->Branch("hepmcp_BCID", _hepmcp_BCID, "hepmcp_BCID[nHepmcp]/I");
  tree->Branch("hepmcp_status", _hepmcp_status, "hepmcp_status[nHepmcp]/I");
  tree->Branch("hepmcp_PDG", _hepmcp_PDG, "hepmcp_PDG[nHepmcp]/I");
  tree->Branch("hepmcp_E", _hepmcp_E, "hepmcp_E[nHepmcp]/F");
  tree->Branch("hepmcp_px", _hepmcp_px, "hepmcp_px[nHepmcp]/F");
  tree->Branch("hepmcp_py", _hepmcp_py, "hepmcp_py[nHepmcp]/F");
  tree->Branch("hepmcp_pz", _hepmcp_pz, "hepmcp_pz[nHepmcp]/F");
  tree->Branch("hepmcp_m1", _hepmcp_m1, "hepmcp_m1[nHepmcp]/I");
  tree->Branch("hepmcp_m2", _hepmcp_m2, "hepmcp_m2[nHepmcp]/I");

  tree->Branch("hepmcp_x1", &_hepmcp_x1, "hepmcp_x1/F");
  tree->Branch("hepmcp_x2", &_hepmcp_x2, "hepmcp_x2/F");
  tree->Branch("hepmcp_Q2", &_hepmcp_Q2, "hepmcp_Q2/F");

   tree->Branch("nTracks", &_nTracks, "nTracks/I");
   tree->Branch("track_ID", _track_ID, "track_ID[nTracks]F");
   tree->Branch("track_trueID", _track_trueID, "track_trueID[nTracks]F");
   tree->Branch("track_px", _track_px, "track_px[nTracks]F");
   tree->Branch("track_py", _track_py, "track_py[nTracks]F");
   tree->Branch("track_pz", _track_pz, "track_pz[nTracks]F");
   tree->Branch("track_dca", _track_dca, "track_dca[nTracks]F");
   tree->Branch("track_dca_2d", _track_dca_2d, "track_dca_2d[nTracks]F");
   tree->Branch("tracks_charge", _track_charge, "tracks_charge[nTracks]/S");
   tree->Branch("track_source", _track_source, "track_source[nTracks]F");
    
   tree->Branch("maxNProjections", &_maxNProjections, "maxNProjections/I");
   tree->Branch("track_ProjTrackID", _track_ProjTrackID, "track_ProjTrackID[maxNProjections]/F");
   tree->Branch("track_ProjLayer", _track_ProjLayer, "track_ProjLayer[maxNProjections]/F");
   tree->Branch("track_TLP_x", _track_TLP_x, "track_TLP_x[maxNProjections]/F");
   tree->Branch("track_TLP_y", _track_TLP_y, "track_TLP_y[maxNProjections]/F");
   tree->Branch("track_TLP_z", _track_TLP_z, "track_TLP_z[maxNProjections]/F");
   tree->Branch("rack_TLP_t", _track_TLP_t, "rack_TLP_t[maxNProjections]/F");
   tree->Branch("track_TLP_true_x", _track_TLP_true_x, "track_TLP_true_x[maxNProjections]/F");
   tree->Branch("track_TLP_true_y", _track_TLP_true_y, "track_TLP_true_y[maxNProjections]/F");
   tree->Branch("track_TLP_true_z", _track_TLP_true_z, "track_TLP_true_z[maxNProjections]/F");

    tree->Branch("B0hits", &_B0hits, "B0hits/I");
    tree->Branch("B0x", _B0x, "B0x[B0hits]/F");
    tree->Branch("B0y", _B0y, "B0y[B0hits]/F");
    tree->Branch("B0z", _B0z, "B0z[B0hits]/F");
    tree->Branch("B0ind", _B0ind, "B0ind[B0hits]/I");
    tree->Branch("B0px", _B0px, "B0px[B0hits]/F");
    tree->Branch("B0py", _B0py, "B0py[B0hits]/F");
    tree->Branch("B0pz", _B0pz, "B0pz[B0hits]/F");
    tree->Branch("B0trPx", _B0trPx, "B0trPx[B0hits]/F");
    tree->Branch("B0trPy", _B0trPy, "B0trPy[B0hits]/F");
    tree->Branch("B0trPz", _B0trPz, "B0trPz[B0hits]/F");
    tree->Branch("B0id", _B0id, "B0id[B0hits]/I");

    tree->Branch("RP1",&_RP1,"RP1/I");
    tree->Branch("RP2",&_RP2,"RP2/I");
    tree->Branch("RPhits",&_RPhits,"RPhits/I");
    tree->Branch("RPx",_RPx,"RPx[RPhits]/F");
    tree->Branch("RPy",_RPy,"RPy[RPhits]/F");
    tree->Branch("RPz",_RPz,"RPz[RPhits]/F");
    tree->Branch("RPind",_RPind,"RPind[RPhits]/I");
    tree->Branch("RPtrPx",_RPtrPx,"RPtrPx[RPhits]/F");
    tree->Branch("RPtrPy",_RPtrPy,"RPtrPy[RPhits]/F");
    tree->Branch("RPtrPz",_RPtrPz,"RPtrPz[RPhits]/F");
    tree->Branch("RPid",_RPid,"RPid[RPhits]/I");
    tree->Branch("RPpx",_RPpx,"RPpx[RPhits]/F");
    tree->Branch("RPpy",_RPpy,"RPpy[RPhits]/F");
    tree->Branch("RPpz",_RPpz,"RPpz[RPhits]/F");
    tree->Branch("RPpid",_RPpid,"RPpid[RPhits]/I");
   
    tree->Branch("ZDChits", &_ZDChits, "ZDChits/I");
    tree->Branch("ZDCx", _ZDCx, "ZDCx[ZDChits]/F");
    tree->Branch("ZDCy", _ZDCy, "ZDCy[ZDChits]/F");
    tree->Branch("ZDCz", _ZDCz, "ZDCz[ZDChits]/F");
    tree->Branch("ZDCpx", _ZDCpx, "ZDCpx[ZDChits]/F");
    tree->Branch("ZDCpy", _ZDCpy, "ZDCpy[ZDChits]/F");
    tree->Branch("ZDCpz", _ZDCpz, "ZDCpz[ZDChits]/F");
    tree->Branch("ZDCtrPx", _ZDCtrPx, "ZDCtrPx[ZDChits]/F");
    tree->Branch("ZDCtrPy", _ZDCtrPy, "ZDCtrPy[ZDChits]/F");
    tree->Branch("ZDCtrPz", _ZDCtrPz, "ZDCtrPz[ZDChits]/F");
    tree->Branch("ZDCid", _ZDCid, "ZDCid[ZDChits]/I");
/*   
    tree->Branch("nMCPart", &_nMCPart, "nMCPart/I");
    tree->Branch("mcpart_ID", _mcpart_ID, "mcpart_ID[nMCPart]/I");
    tree->Branch("mcpart_ID_parent", _mcpart_ID_parent, "mcpart_ID_parent[nMCPart]/I");
    tree->Branch("mcpart_PDG", _mcpart_PDG, "mcpart_PDG[nMCPart]/I");
    tree->Branch("mcpart_E", _mcpart_E, "mcpart_E[nMCPart]/F");
    tree->Branch("mcpart_px", _mcpart_px, "mcpart_px[nMCPart]/F");
    tree->Branch("mcpart_py", _mcpart_py, "mcpart_py[nMCPart]/F");
    tree->Branch("mcpart_pz", _mcpart_pz, "mcpart_pz[nMCPart]/F");
    tree->Branch("mcpart_BCID", _mcpart_BCID, "mcpart_BCID[nMCPart]/I");
*/
    tree->Branch("tower_FHCAL_N", &_nTowers_FHCAL, "tower_FHCAL_N/I");
    tree->Branch("tower_FHCAL_E", _tower_FHCAL_E, "tower_FHCAL_E[tower_FHCAL_N]/F");
    tree->Branch("tower_FHCAL_iEta", _tower_FHCAL_iEta, "tower_FHCAL_iEta[tower_FHCAL_N]/I");
    tree->Branch("tower_FHCAL_iPhi", _tower_FHCAL_iPhi, "tower_FHCAL_iPhi[tower_FHCAL_N]/I");
    tree->Branch("tower_FHCAL_trueID", _tower_FHCAL_trueID, "tower_FHCAL_trueID[tower_FHCAL_N]/I");

    tree->Branch("tower_BECAL_N", &_nTowers_BECAL, "tower_BECAL_N/I");
    tree->Branch("tower_BECAL_E", _tower_BECAL_E, "tower_BECAL_E[tower_BECAL_N]/F");
    tree->Branch("tower_BECAL_iEta", _tower_BECAL_iEta, "tower_BECAL_iEta[tower_BECAL_N]/I");
    tree->Branch("tower_BECAL_iPhi", _tower_BECAL_iPhi, "tower_BECAL_iPhi[tower_BECAL_N]/I");
    tree->Branch("tower_BECAL_trueID", _tower_BECAL_trueID, "tower_BECAL_trueID[tower_BECAL_N]/I");

    tree->Branch("tower_HCALIN_N", &_nTowers_HCALIN, "tower_HCALIN_N/I");
    tree->Branch("tower_HCALIN_E", _tower_HCALIN_E, "tower_HCALIN_E[tower_HCALIN_N]/F");
    tree->Branch("tower_HCALIN_iEta", _tower_HCALIN_iEta, "tower_HCALIN_iEta[tower_HCALIN_N]/I");
    tree->Branch("tower_HCALIN_iPhi", _tower_HCALIN_iPhi, "tower_HCALIN_iPhi[tower_HCALIN_N]/I");
    tree->Branch("tower_HCALIN_trueID", _tower_HCALIN_trueID, "tower_HCALIN_trueID[tower_HCALIN_N]/I"); 
 
    tree->Branch("tower_HCALOUT_N", &_nTowers_HCALOUT, "tower_HCALOUT_N/I");
    tree->Branch("tower_HCALOUT_E", _tower_HCALOUT_E, "tower_HCALOUT_E[tower_HCALOUT_N]/F");
    tree->Branch("tower_HCALOUT_iEta", _tower_HCALOUT_iEta, "tower_HCALOUT_iEta[tower_HCALOUT_N]/I");
    tree->Branch("tower_HCALOUT_iPhi", _tower_HCALOUT_iPhi, "tower_HCALOUT_iPhi[tower_HCALOUT_N]/I");
    tree->Branch("tower_HCALOUT_trueID", _tower_HCALOUT_trueID, "tower_HCALOUT_trueID[tower_HCALOUT_N]/I");

    tree->Branch("tower_EHCAL_N", &_nTowers_EHCAL, "tower_EHCAL_N/I");
    tree->Branch("tower_EHCAL_E", _tower_EHCAL_E, "tower_EHCAL_E[tower_EHCAL_N]/F");
    tree->Branch("tower_EHCAL_iEta", _tower_EHCAL_iEta, "tower_EHCAL_iEta[tower_EHCAL_N]/I");
    tree->Branch("tower_EHCAL_iPhi", _tower_EHCAL_iPhi, "tower_EHCAL_iPhi[tower_EHCAL_N]/I");
    tree->Branch("tower_EHCAL_trueID", _tower_EHCAL_trueID, "tower_EHCAL_trueID[tower_EHCAL_N]/I");
   
    tree->Branch("tower_FEMC_N", &_nTowers_FEMC, "tower_FEMC_N/I");
    tree->Branch("tower_FEMC_E", _tower_FEMC_E, "tower_FEMC_E[tower_FEMC_N]/F");
    tree->Branch("tower_FEMC_iEta", _tower_FEMC_iEta, "tower_FEMC_iEta[tower_FEMC_N]/I");
    tree->Branch("tower_FEMC_iPhi", _tower_FEMC_iPhi, "tower_FEMC_iPhi[tower_FEMC_N]/I");
    tree->Branch("tower_FEMC_trueID", _tower_FEMC_trueID, "tower_FEMC_trueID[tower_FEMC_N]/I");

    tree->Branch("tower_EEMC_N", &_nTowers_EEMC, "tower_EEMC_N/I");
    tree->Branch("tower_EEMC_E", _tower_EEMC_E, "tower_EEMC_E[tower_EEMC_N]/F");
    tree->Branch("tower_EEMC_iEta", _tower_EEMC_iEta, "tower_EEMC_iEta[tower_EEMC_N]/I");
    tree->Branch("tower_EEMC_iPhi", _tower_EEMC_iPhi, "tower_EEMC_iPhi[tower_EEMC_N]/I");
    tree->Branch("tower_EEMC_trueID", _tower_EEMC_trueID, "tower_EEMC_trueID[tower_EEMC_N]/I");

/*    tree->Branch("vertex_true_x", &_vertex_true_x, "vertex_true_x/F");
    tree->Branch("vertex_true_y", &_vertex_true_y, "vertex_true_y/F");
    tree->Branch("vertex_true_z", &_vertex_true_z, "vertex_true_z/F");
  */
}
//**************************************88
void excl_ana::initializeVariables()
{

  _EMJpsi_px = _EMJpsi_py = _EMJpsi_pz = _EMJpsi_e = -1000.;
  _EPJpsi_px = _EPJpsi_py = _EPJpsi_pz = _EPJpsi_e = -1000.;
  _E_px = _E_py = _E_pz = _E_e = -1000.;
  _P_px = _P_py = _P_pz = _P_e = -1000.;

  for (int k=0; k<_maxNHepmcp; k++){
	
	_hepmcp_BCID[k] = _hepmcp_status[k] = _hepmcp_PDG[k] = -1000.;
        _hepmcp_E[k] = _hepmcp_px[k] = _hepmcp_py[k] = _hepmcp_pz[k] = -1000.;
	_hepmcp_m1[k] = _hepmcp_m2[k] =  -1000. ; 
  }

  for (int k=0; k<_maxNTracks; k++){
    _track_ID[k] = _track_trueID[k] = _track_px[k] = _track_py[k] = _track_pz[k] =  -1000.;
    _track_dca[k] = _track_dca_2d[k]  =  -1000.;
  }
   for (int k=0; k<_maxNProjections; k++){
	_track_ProjTrackID[k] = _track_ProjLayer[k] = _track_TLP_x[k] = _track_TLP_y[k] = _track_TLP_z[k] =  -1000.;
	_track_TLP_t[k] = _track_TLP_true_x[k] = _track_TLP_true_y[k] =  _track_TLP_true_z[k] = -1000.;
   }
}


 int excl_ana::GetProjectionIndex(std::string projname)
 {
    if (projname.find("HCALIN") != std::string::npos)
      return 1;
    else if (projname.find("HCALOUT") != std::string::npos)
      return 2;
    else if (projname.find("CEMC") != std::string::npos)
      return 3;
    else
      return -1;
   return -1;
 }

std::string excl_ana::GetProjectionNameFromIndex(int projindex)
 {
    switch (projindex)
   {
      case 1:
       return "HCALIN";
      case 2:
        return "HCALOUT";
     case 3:
        return "CEMC";
      default:
        return "NOTHING";
    }
 }
//=================================================================
// B0 Hits

int excl_ana::process_B0(PHCompositeNode* topNode)
{
  ostringstream nodename;
  nodename.str("");
  nodename << "G4HIT_" << "b0Truth";
  
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  int layer=-1;
  _B0hits = 0;
  if (hits) {

    PHG4HitContainer::ConstRange hit_range = hits->getHits();


    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {
	
	
        PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
	PHG4Particle* g4particle = truthinfo->GetParticle(hit_iter->second->get_trkid());
	if(!g4particle)	continue;		
	
     	_B0x[_B0hits] = hit_iter->second->get_x(0);
	_B0y[_B0hits] = hit_iter->second->get_y(0);
	_B0z[_B0hits] = hit_iter->second->get_z(0);

	_B0px[_B0hits] = hit_iter->second->get_px(0);
        _B0py[_B0hits] = hit_iter->second->get_py(0);
        _B0pz[_B0hits] = hit_iter->second->get_pz(0);
 	_B0id[_B0hits]  = hit_iter->second->get_hit_id();

	_B0pid[_B0hits]=g4particle->get_pid();
	
	if(TMath::Abs(_B0z[_B0hits] - 592.0)<5) {layer = 1;}
        if(TMath::Abs(_B0z[_B0hits] - 616.0)<5) {layer = 2;}
        if(TMath::Abs(_B0z[_B0hits] - 640.0)<5) {layer = 3;}
        if(TMath::Abs(_B0z[_B0hits] - 663.0)<5) {layer = 4;}
	
	_B0ind[_B0hits] = layer;
  	PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();  

        for (PHG4TruthInfoContainer::ConstIterator iter = range.first;  iter != range.second; ++iter){
                const PHG4Particle *truth = iter->second;
                if(truth->get_pid() == 2212){   
                    _B0trPx[_B0hits] = truth->get_px();
                    _B0trPy[_B0hits] = truth->get_py();
                    _B0trPz[_B0hits] = truth->get_pz();
                }               
         } // end PHG4TruthInfoContainer iterator
	
    } _B0hits++;
	
 }
  return Fun4AllReturnCodes::EVENT_OK;
}


//=================================================================
//Roman Pots

int excl_ana::process_RomanPots(PHCompositeNode* topNode, int rPot){
  
  ostringstream nodename;
  nodename.str("");
  nodename << "G4HIT_" << "rpTruth";
  
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
   _RPhits = 0;
  if (hits) {
 
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    int counter = 0;
    int counter2 = 0;

     for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {

	_RPx[_RPhits] = hit_iter->second->get_x(0);
        _RPy[_RPhits] = hit_iter->second->get_y(0);
        _RPz[_RPhits] = hit_iter->second->get_z(0);
	_RPpz[_RPhits] = hit_iter->second->get_pz(0);
	_RPpx[_RPhits] = hit_iter->second->get_px(0);
	_RPpy[_RPhits] = hit_iter->second->get_py(0);
	_RPid[_RPhits] = hit_iter->second->get_hit_id();
	if(TMath::Abs(_RPz[_RPhits] - 2600.0)<50) {rPot = 1; counter++;}
	if(TMath::Abs(_RPz[_RPhits] - 2800.0)<50) {rPot = 2; counter2++;}
	_RPind[_RPhits ] = rPot;
  
	PHG4Particle* g4particle = truthinfo->GetParticle(hit_iter->second->get_trkid());
        _RPpid[_RPhits] = 0;
	_RPpid[_RPhits] = g4particle->get_pid();
  
	for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter){

	  const PHG4Particle *truth = iter->second;
	  if(truth->get_pid() == 2212){
	    _RPtrPx[_RPhits] = truth->get_px();
            _RPtrPy[_RPhits] = truth->get_py();
            _RPtrPz[_RPhits] = truth->get_pz();
	  }
	}
	_RPhits++; 
      }
      _RP1 = counter;
      _RP2 = counter2;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}



int excl_ana::process_ZDC(PHCompositeNode* topNode){

  _ZDChits = 0;

  ostringstream nodename;
  nodename.str("");
  nodename << "G4HIT_" << "ZDC";
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

  if (hits) {
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
   
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
 
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {	

	_ZDCx[_ZDChits] = (Float_t)hit_iter->second->get_x(0);
	_ZDCy[_ZDChits] = (Float_t)hit_iter->second->get_y(0);
	_ZDCz[_ZDChits] = (Float_t)hit_iter->second->get_z(0);

	_ZDCpz[_ZDChits] = hit_iter->second->get_pz(0);
        _ZDCpx[_ZDChits] = hit_iter->second->get_px(0);
        _ZDCpy[_ZDChits] = hit_iter->second->get_py(0);
        _ZDCid[_ZDChits] = hit_iter->second->get_hit_id();

	PHG4Particle* g4particle = truthinfo->GetParticle(hit_iter->second->get_trkid());

	_ZDCpid[_ZDChits]=0;
	_ZDCpid[_ZDChits]=g4particle->get_pid();

	
	 for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter){
	   const PHG4Particle *truth = iter->second;
	   if(truth->get_pid() == 2212){
	     _ZDCtrPx[_ZDChits] = truth->get_px();
             _ZDCtrPy[_ZDChits] = truth->get_py();
             _ZDCtrPz[_ZDChits] = truth->get_pz();
	   } 	   
	 }//end truth container

	_ZDChits++; 
    }//end for hit container
	
  }//end hits
  return Fun4AllReturnCodes::EVENT_OK;
} // end szdc







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



using namespace std;


//____________________________________________________________________________..
//excl_ana::excl_ana(const std::string &name):
// SubsysReco(name)
//{
//  std::cout << "excl_ana::excl_ana(const std::string &name) Calling ctor" << std::endl;
//}


excl_ana::excl_ana(const std::string &name):
 SubsysReco(name)
{

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

  set_do_PROJECTIONS(false);
  initializeVariables();
  initializeTrees();


}


void excl_ana::SetOutputFile(TString outfilename) {

        mOutputFileName = outfilename;

}


//____________________________________________________________________________..
excl_ana::~excl_ana()
{
  delete tree;
if(_caloevalstackCEMC) delete _caloevalstackCEMC;
if(_caloevalstackFEMC) delete _caloevalstackFEMC;
if(_caloevalstackEEMC) delete _caloevalstackEEMC;
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


PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");


PHG4Particle* primary;

int pid, barcode;
//cout << "======================================================" << endl;
for(int kk=0;kk<truthinfocontainer->GetNumPrimaryVertexParticles();kk++){
	primary = truthinfocontainer->GetPrimaryParticle(kk+1);
	pid     = primary->get_pid();
	barcode = primary->get_barcode();

	if(pid == 11 && barcode == 10009){
		_EMJpsi_px = primary->get_px();
		_EMJpsi_py = primary->get_py();
		_EMJpsi_pz = primary->get_pz();
		_EMJpsi_e  = primary->get_e();
	}
	if(pid == -11 && barcode == 10008){
		_EPJpsi_px = primary->get_px();
                _EPJpsi_py = primary->get_py();
                _EPJpsi_pz = primary->get_pz();
                _EPJpsi_e  = primary->get_e();
	}
	if(pid == 2212  && barcode == 10007){
                _P_px = primary->get_px();
                _P_py = primary->get_py();
                _P_pz = primary->get_pz();
                _P_e  = primary->get_e();
        }

	if(pid == 11 && barcode == 10002){
                _E_px = primary->get_px();
                _E_py = primary->get_py();
                _E_pz = primary->get_pz();
                _E_e  = primary->get_e();
        }	
}

_nMCPart = 0;

if (truthinfocontainer){
   PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetParticleRange();
  
    for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr){
     PHG4Particle* g4particle = truth_itr->second;
     if (!g4particle) continue;
     if (g4particle->get_e() < 1.) continue;    
 
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
     return 0;
}



// ========= lopping over HEPMC Particles =============================================//

PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

if (hepmceventmap){

  int _nHepmcp = 0 ;
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
int _nTracks = 0;
int _nProjections = 0;  
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
  tree->Branch("EMJpsi_px", &_EMJpsi_px, "EMJpsi_px/F");
  tree->Branch("EMJpsi_py", &_EMJpsi_py, "EMJpsi_py/F");
  tree->Branch("EMJpsi_pz", &_EMJpsi_pz, "EMJpsi_pz/F");
  tree->Branch("EMJpsi_e",  &_EMJpsi_e,  "EMJpsi_e/F");

  tree->Branch("EPJpsi_px", &_EPJpsi_px, "EPJpsi_px/F");
  tree->Branch("EPJpsi_py", &_EPJpsi_py, "EPJpsi_py/F");
  tree->Branch("EPJpsi_pz", &_EPJpsi_pz, "EPJpsi_pz/F");
  tree->Branch("EPJpsi_e",  &_EPJpsi_e,  "EPJpsi_e/F");

  tree->Branch("E_px", &_E_px, "E_px/F");
  tree->Branch("E_py", &_E_py, "E_py/F");
  tree->Branch("E_pz", &_E_pz, "E_pz/F");
  tree->Branch("E_e",  &_E_e,  "E_e/F");

  tree->Branch("P_px", &_P_px, "P_px/F");
  tree->Branch("P_py", &_P_py, "P_py/F");
  tree->Branch("P_pz", &_P_pz, "P_pz/F");
  tree->Branch("P_e",  &_P_e,  "P_e/F");

  tree->Branch("maxNHepmcp", &_maxNHepmcp, "maxNHepmcp/I");
  tree->Branch("hepmcp_BCID", _hepmcp_BCID, "hepmcp_BCID[maxNHepmcp]/I");
  tree->Branch("hepmcp_status", _hepmcp_status, "hepmcp_status[maxNHepmcp]/I");
  tree->Branch("hepmcp_PDG", _hepmcp_PDG, "hepmcp_PDG[maxNHepmcp]/I");
  tree->Branch("hepmcp_E", _hepmcp_E, "hepmcp_E[maxNHepmcp]/F");
  tree->Branch("hepmcp_px", _hepmcp_px, "hepmcp_px[maxNHepmcp]/F");
  tree->Branch("hepmcp_py", _hepmcp_py, "hepmcp_py[maxNHepmcp]/F");
  tree->Branch("hepmcp_pz", _hepmcp_pz, "hepmcp_pz[maxNHepmcp]/F");
  tree->Branch("hepmcp_m1", _hepmcp_m1, "hepmcp_m1[maxNHepmcp]/I");
  tree->Branch("hepmcp_m2", _hepmcp_m2, "hepmcp_m2[maxNHepmcp]/I");

  tree->Branch("hepmcp_x1", &_hepmcp_x1, "hepmcp_x1/F");
  tree->Branch("hepmcp_x2", &_hepmcp_x2, "hepmcp_x2/F");
  tree->Branch("hepmcp_Q2", &_hepmcp_Q2, "hepmcp_Q2/F");

   tree->Branch("nTracks", &_maxNTracks, "nTracks/I");
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
   
    tree->Branch("nMCPart", &_nMCPart, "nMCPart/I");
    tree->Branch("mcpart_ID", _mcpart_ID, "mcpart_ID[nMCPart]/I");
    tree->Branch("mcpart_ID_parent", _mcpart_ID_parent, "mcpart_ID_parent[nMCPart]/I");
    tree->Branch("mcpart_PDG", _mcpart_PDG, "mcpart_PDG[nMCPart]/I");
    tree->Branch("mcpart_E", _mcpart_E, "mcpart_E[nMCPart]/F");
    tree->Branch("mcpart_px", _mcpart_px, "mcpart_px[nMCPart]/F");
    tree->Branch("mcpart_py", _mcpart_py, "mcpart_py[nMCPart]/F");
    tree->Branch("mcpart_pz", _mcpart_pz, "mcpart_pz[nMCPart]/F");
    tree->Branch("mcpart_BCID", _mcpart_BCID, "mcpart_BCID[nMCPart]/I");

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

  Int_t layer=-1;
  _B0hits = 0;
  if (hits) {
  PHG4HitContainer::ConstRange hit_range = hits->getHits();

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();  

  for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {

     	_B0x[_B0hits] = (Float_t)hit_iter->second->get_x(0);
        _B0y[_B0hits] = (Float_t)hit_iter->second->get_y(0);
        _B0z[_B0hits] = (Float_t)hit_iter->second->get_z(0);

        _B0pz[_B0hits] = hit_iter->second->get_pz(0);
        _B0px[_B0hits] = hit_iter->second->get_px(0);
        _B0py[_B0hits] = hit_iter->second->get_py(0);
        _B0id[_B0hits] = hit_iter->second->get_hit_id();

	PHG4Particle* g4particle = truthinfo->GetParticle(hit_iter->second->get_trkid());
	_B0pid[_B0hits]=0;
	_B0pid[_B0hits]=g4particle->get_pid();

	if(TMath::Abs(_B0z[_B0hits] - 592.0)<5) {layer = 1;}
	if(TMath::Abs(_B0z[_B0hits] - 616.0)<5) {layer = 2;}
	if(TMath::Abs(_B0z[_B0hits] - 640.0)<5) {layer = 3;}
	if(TMath::Abs(_B0z[_B0hits] - 663.0)<5) {layer = 4;}
	
	_B0ind[_B0hits] = layer;

	for (PHG4TruthInfoContainer::ConstIterator iter = range.first;  iter != range.second; ++iter){
	
	    	const PHG4Particle *truth = iter->second;
	  	if(truth->get_pid() == 2212){	
		    _B0trPx[_B0hits] = truth->get_px();
		    _B0trPy[_B0hits] = truth->get_py();
		    _B0trPz[_B0hits] = truth->get_pz();
		} 		
	 } // end PHG4TruthInfoContainer iterator
	  _B0hits++;
	}
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







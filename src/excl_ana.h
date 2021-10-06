// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EXCL_ANA_H
#define EXCL_ANA_H

#include <fun4all/SubsysReco.h>

#include <string>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"


class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TTree;
class TNtuple;
class CaloEvalStack;
class CaloRawClusterEval;
class RawClusterContainer;
class SvtxTrackMap;
class SvtxEvalStack;
class SvtxTrackEval;
class PHG4TruthInfoContainer;


class excl_ana : public SubsysReco
{
 public:
 enum class TrackSource_t : unsigned short{
   all = 0,
   inner =1 
 };

CaloEvalStack *_caloevalstack;
CaloEvalStack* _caloevalstackFEMC;
CaloEvalStack* _caloevalstackEEMC;
CaloEvalStack* _caloevalstackCEMC;

SvtxEvalStack *_svtxEvalStack;
  excl_ana(const std::string &name = "excl_ana");
 // excl_ana(const std::string &name = "Diff_Tagg_ana", const std::string &fname = "MyNtuple.root");

  virtual ~excl_ana();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int process_g4hits(PHCompositeNode *);

  int process_RomanPots(PHCompositeNode *, int);
  int process_B0(PHCompositeNode *);
  int process_ZDC(PHCompositeNode *);
 
 int process_ClusterCalo(PHCompositeNode*, std::string);

  void set_do_PROJECTIONS(bool b) { _do_PROJECTIONS = b; } 

  void SetOutputFile(TString);
 private:

  TTree *tree;

  float _EMJpsi_px, _EMJpsi_py, _EMJpsi_pz, _EMJpsi_e;
  float _EPJpsi_px, _EPJpsi_py, _EPJpsi_pz, _EPJpsi_e;
  float _E_px, _E_py, _E_pz, _E_e;
  float _P_px, _P_py, _P_pz, _P_e;

  int _maxNHepmcp = 50; 
  int* _hepmcp_BCID = new int[_maxNHepmcp];
  int*  _hepmcp_status = new int[_maxNHepmcp];
  int* _hepmcp_PDG = new int[_maxNHepmcp];
  float* _hepmcp_E = new float[_maxNHepmcp];
  float* _hepmcp_px = new float[_maxNHepmcp];
  float* _hepmcp_py = new float[_maxNHepmcp];
  float* _hepmcp_pz = new float[_maxNHepmcp];
  int* _hepmcp_m1 = new int[_maxNHepmcp];
  int* _hepmcp_m2 = new int[_maxNHepmcp];
  float _hepmcp_x1;
  float _hepmcp_x2;
  float _hepmcp_Q2;

  int _maxNTracks = 50;

  float* _track_ID = new float[_maxNTracks];
  float* _track_trueID = new float[_maxNTracks];
  float*  _track_px = new float[_maxNTracks];
  float* _track_py = new float[_maxNTracks];
  float* _track_pz = new float[_maxNTracks];
  float* _track_dca = new float[_maxNTracks];
  float* _track_dca_2d = new float[_maxNTracks];
  unsigned short* _track_source = new unsigned short[_maxNTracks];

  int _maxNProjections =  50;
  float* _track_ProjTrackID = new float[_maxNProjections];
  int* _track_ProjLayer = new int[_maxNProjections];
  float* _track_TLP_x = new float[_maxNProjections];
  float* _track_TLP_y = new float[_maxNProjections];
  float* _track_TLP_z = new float[_maxNProjections];
  float* _track_TLP_t = new float[_maxNProjections];
  float* _track_TLP_true_x = new float[_maxNProjections];
  float* _track_TLP_true_y = new float[_maxNProjections];
  float* _track_TLP_true_z = new float[_maxNProjections];
  float* _track_TLP_true_t = new float[_maxNProjections];
  short* _track_charge = new short[_maxNProjections]; 

   int _B0hits = 100;
    
   float* _B0x =  new float[_B0hits];
   float* _B0y =  new float[_B0hits];
   float* _B0z =  new float[_B0hits];
   float* _B0px =  new float[_B0hits];
   float* _B0py =  new float[_B0hits];
   float* _B0pz =   new float[_B0hits];
   float* _B0trPx =  new float[_B0hits];
   float* _B0trPy =  new float[_B0hits];
   float* _B0trPz =  new float[_B0hits];
   int*  _B0pid =  new int[_B0hits];  
   int*  _B0id =  new int[_B0hits];  
   int*  _B0ind  =  new int[_B0hits];

   int _RP1  = 100;
   int _RP2  = 100;
   int _RPhits  = 100;
    
   float* _RPx = new float[_RPhits];
   float* _RPy = new float[_RPhits];
   float* _RPz = new float[_RPhits];
   int*   _RPind = new int[_RPhits];
   float* _RPpx = new float[_RPhits];
   float* _RPpy = new float[_RPhits];
   float* _RPpz = new float[_RPhits];
   float* _RPtrPx = new float[_RPhits];
   float* _RPtrPy = new float[_RPhits];
   float* _RPtrPz = new float[_RPhits];
   int*   _RPid = new int[_RPhits];
   int*  _RPpid = new int[_RPhits];


   int _ZDChits = 100;

   float* _ZDCx    =  new float[_ZDChits];
   float* _ZDCy    =  new float[_ZDChits];
   float* _ZDCz    =  new float[_ZDChits];
   float* _ZDCpx   =  new float[_ZDChits];
   float* _ZDCpy   =  new float[_ZDChits];
   float* _ZDCpz   =  new float[_ZDChits];
   float* _ZDCtrPx =  new float[_ZDChits];
   float* _ZDCtrPy =  new float[_ZDChits];
   float* _ZDCtrPz =  new float[_ZDChits];
   int*   _ZDCpid  =  new int[_ZDChits];
   int*   _ZDCid   =  new int[_ZDChits];

   int _maxNMCPart = 100;
   int _nMCPart    = 100;
   int* _mcpart_ID = new int[_maxNMCPart];
   int* _mcpart_ID_parent = new int[_maxNMCPart];
   int* _mcpart_PDG = new int[_maxNMCPart];
   float* _mcpart_E = new float[_maxNMCPart];
   float* _mcpart_px = new float[_maxNMCPart];
   float* _mcpart_py = new float[_maxNMCPart];
   float* _mcpart_pz = new float[_maxNMCPart];
   int*_mcpart_BCID = new int[_maxNMCPart];

   TString mOutputFileName;


   float _vertex_true_x;
   float _vertex_true_y;
    float _vertex_true_z;

  void initializeVariables();
  void initializeTrees();
  bool _do_PROJECTIONS;
  std::string GetProjectionNameFromIndex(int projindex);  ///< return track projection layer name from projection index (see GetProjectionIndex)
    int GetProjectionIndex(std::string projname);           ///< return track projection index for given track projection layer
 
protected:

  std::string detector;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  
  unsigned long long int event_itt;
  gsl_rng* m_RandomGenerator;


};

#endif // DIFF_TAGG_ANA_H

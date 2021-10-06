#ifndef MACRO_FUN4ALLG4RUNEVALUATORS_C
#define MACRO_FUN4ALLG4RUNEVALUATORS_C

#include <exclusive/excl_ana.h>
#include <dirent.h>
#include <fstream>
#include <map>
#include <stdlib.h>

#include <GlobalVariables.C>

#include <G4Setup_EICDetector.C>
#include <G4_EventEvaluator.C>
#include <G4_FwdJets.C>
#include <G4_Global.C>
#include <G4_Input.C>
#include <G4_Production.C>
#include <G4_User.C>

#include <fun4all/Fun4AllServer.h>

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(excl_ana.so)


int Fun4All_exclusive(
    const int nEvents,
    const string &inputFile,
    const string &outputFile)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  	Input::READHITS = true;
  	// If you are specifying a single DST
	INPUTREADHITS::filename[0] = inputFile;
  	// If you are specifying a filelist
//	INPUTREADHITS::listfile[0] = inputFile;

  InputInit();

  excl_ana* excRecon = new excl_ana();
  excRecon->SetOutputFile(outputFile);
  se->registerSubsystem(excRecon);

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  //-----
  // Run
  //-----

  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  cout << "All done" << endl;
  delete se;

  gSystem->Exit(0);
  return 0;
}

#endif

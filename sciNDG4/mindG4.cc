/// ---------------------------------------------------------------------------
///  \file   mindG4.cc 
///  \brief  Main program.
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>    
///  \date     16 Jun 2008
///  \version  $Id: mindG4.cc 543 2014-11-01 23:03:10Z  $
///
///  Copyright (c) 2008, 2009 -- IFIC Neutrino Group
/// ---------------------------------------------------------------------------

#include "MindConfigService.h"
#include "MindParamStore.h"
#include "MindDetectorConstruction.h"
#include "MindPrimaryGeneration.h"
#include "MindRunAction.h"
#include "MindEventAction.h"
#include "MindTrackingAction.h"
#include "MindSteppingAction.h"
// #include "MindQGSP_BERT.h"
#include "MindUtils.h"

#include <G4RunManager.hh>
#include <G4UImanager.hh>   
#include <G4UIterminal.hh>  
#include <G4UItcsh.hh> 
#include <G4VisExecutive.hh>
#include <G4PhysListFactory.hh>
#include <globals.hh>
#include <map>

#include <time.h>

void PrintUsage()
{
  G4cerr << "\nUsage: mindG4 [mode] [mode_options] <config_file>\n" << G4endl;
  G4cerr << "Available modes:" << G4endl;
  G4cerr << "   [-batch(b)]               : "
	 << "batch mode" << G4endl;
  G4cerr << "    -visual(v) [<vis_macro>] : "
	 << "visual mode (default macro: vis.mac)" << G4endl;
  G4cerr << G4endl;
  exit(-1);
}



int main(int argc, char** argv) 
{
  // command-line options ...........................................
  
  G4bool   visual     = false;
  G4String vis_macro  = "vis.mac"; // default visual macro
  G4String config_file;  
  
  if (argc < 2 || argc > 4) {  // too many or too few arguments
    PrintUsage();
  }
  else if (argc == 2) {        // (implicit) batch mode 
    config_file = argv[1];
  }
  else if (argc > 2) {         // explicit mode (either batch or visual)
    G4String mode = argv[1];
    config_file   = argv[2];
    if (mode == "-visual" || mode == "-v") {  // visual mode
      visual = true;
      if (argc == 4) vis_macro = argv[3];
    }
  }

  // Initialization of required classes .............................
  
  // run-time configuration service
  MindConfigService::Instance().SetConfigFile(config_file);
  
  // Geant4 run manager
  G4RunManager* runMgr = new G4RunManager();
  
  // geometry, generation and physics
  runMgr->SetUserInitialization(new MindDetectorConstruction);
  // runMgr->SetUserInitialization(new MindQGSP_BERT);//QGSP_BERT);
  G4PhysListFactory factory;
  G4VModularPhysicsList* physlist = factory.GetReferencePhysList("QGSP_BERT_EMV");
  runMgr->SetUserInitialization(physlist);
  runMgr->SetUserAction(new MindPrimaryGeneration);
  
  // optional user actions
  runMgr->SetUserAction(new MindRunAction);
  runMgr->SetUserAction(new MindEventAction);
  runMgr->SetUserAction(new MindTrackingAction);
  runMgr->SetUserAction(new MindSteppingAction);

  //Set seed and save it in run.
  G4int seed = MindConfigService::Instance().Job().GetIParam("random_seed");
  if ( seed < 0 )
    CLHEP::HepRandom::setTheSeed( time(0) );
  else
    CLHEP::HepRandom::setTheSeed( seed );
  
  //set bhep output.
  MindUtils::OpenOutputDst( MindConfigService::Instance().Job().GetSParam("output_dst"),
			    seed );

  // VISUAL MODE ....................................................
  if (visual) {
    
    // initialize visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // initialize G4 kernel
    runMgr->Initialize();
    
    // get the pointer to the UI manager and set verbosities
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");
    
    // G4UIterminal is a (dumb) terminal.
    G4UIsession *session = 0;
    session = new G4UIterminal(new G4UItcsh);
    
    // execute visualization macro and start session
    UI->ApplyCommand("/control/execute " + vis_macro);
    session->SessionStart();
    
    delete session;
    delete visManager;
  }
  
  // BATCH MODE .....................................................
  else {
    
    // initialize G4 kernel
    runMgr->Initialize();
    
    // get the pointer to the UI manager and set verbosities
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");
    
    // start the run
    G4int numEvents = 
      MindConfigService::Instance().Job().GetIParam("number_events");
    runMgr->BeamOn(numEvents);
  }

  // close bhep writer
  MindUtils::CloseOutputDst();
  
  delete runMgr;
  return 0;
} 

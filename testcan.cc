// $Id: task2.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file
 * @brief Main program.
 */

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4Version.hh"

#include "G4VisExecutive.hh"
#if  G4VERSION_NUMBER>=930
#include "G4UIExecutive.hh"
#else
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#endif

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PhysicsList.hh"

#include "EventAction.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"

#include "HistoManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <math.h>
#include "QGSP_BERT_HP.hh"
#include "G4VUserPhysicsList.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


/*!
\brief Main program

\callgraph

*/

//class SteppingAction;

int main(int argc,char** argv)
{
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4double seed = time(NULL);
    CLHEP::HepRandom::setTheSeed(seed);
    
    //G4double CrossArray[20];
    
  // Run manager
  G4RunManager * runManager = new G4RunManager();

  // mandatory Initialization classes 
  //G4VUserDetectorConstruction* detector = new DetectorConstruction();
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  G4VUserPhysicsList* physics = new PhysicsList();
  runManager->SetUserInitialization(physics);
    //runManager->SetUserInitialization(new PhysicsList);
  //G4VModularPhysicsList* physicsList = new QGSP_BERT_HP;

  //G4VUserPhysicsList* physicsList = new QGSP_BERT_HP;

  // physicsList->SetVerboseLevel(0);
  // runManager->SetUserInitialization(physicsList);
   
  // mandatory User Action classes
  PrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);
    
  HistoManager * histoMan = new HistoManager();


  //Optional User Action classes
  //Event action (handles for beginning / end of event)
  EventAction* event_action = new EventAction(histoMan);
  runManager->SetUserAction( event_action );
    
    SteppingAction * step_action = new SteppingAction(event_action);
  runManager->SetUserAction( step_action );
    
  RunAction * run_action = new RunAction(histoMan,gen_action,detector);
    runManager->SetUserAction(run_action);
    


  // Initialize G4 kernel

    //     runManager->Initialize();
      
  //Initilize the visualization manager
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();
     
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if (argc!=1) {  // batch mode  
	  //command line contains name of the macro to execute
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
      //UImanager->ApplyCommand("/testcan/phys/SelectPhysics QGSP_BERT_HP");
      //runManager->Initialize();
 }
  else {           // interactive mode : define UI session
     
#if  G4VERSION_NUMBER>=930
	  //New since G4 9.3: UI executive setups up
	  //correct UI depending on env variables
	  G4UIExecutive * ui = new G4UIExecutive(argc,argv);
	  //If UI has graphics execute special macro: opens OpenGL Qt driver
	//  if (ui->IsGUI())
	//	  UImanager->ApplyCommand("/control/execute visQt.mac");
//	  else
		  UImanager->ApplyCommand("/control/execute vis.mac");
#else
	  //Older versions of G4: UI selected by user
  #ifdef G4UI_USE_TCSH
	  G4UIsession * ui = new G4UIterminal(new G4UItcsh);
  #else
	  G4UIsession * ui = new G4UIterminal();
  #endif
	  UImanager->ApplyCommand("/control/execute vis.mac");
#endif
	  ui->SessionStart();
	  delete ui;
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
    
    //G4cout << EventAction::GetnEvents() << " events processes \n" ;
    //G4cout << EventAction::GetnScatteringEvents() << " events interacted \n" ;
    //G4cout << EventAction::GetnEdepEvents() << " events deposited energy \n" ;
    
   G4double scat = EventAction::GetnScatteringEvents();
    G4double inc = EventAction::GetnEvents();
    //G4double edepevents = EventAction::GetnEdepEvents();
  
    G4double avo = 6.022e+23;
    G4double ratio;
    ratio = (scat/inc);
    //G4double csection = (ratio*1e+24)/avo;
    
    
    G4double nlog = double (log(1-ratio));
    G4double crossSect = -nlog/avo;
    
    //G4cout << "\n\n" << "Cross Section: " << crossSect << " cm2" << G4endl;
            //  G4cout << "               " << crossSect*1e+24 << " barns" << G4endl;
    

				// Using Arrays to store cross section calculations

    /* G4double ener = 0.5;
     for(int i=0; i<20 ; i++){
        
      G4cout << run_action->GetCrossArray(i) << " +/- " << run_action->GeterrArray(i) << " barns  ||  " << run_action->GetHitArray(i) << " +/- " << sqrt(run_action->GetHitArray(i)) << " hits  ||  " << ener << " MeV" << G4endl;
      ener += 0.5;
     }*/

				// Using vectors to store cross section calculations

	//G4cout << "number of runs: " << run_action->GetNumberOfRuns() << G4endl;
	//G4cout << "size of cross section vector: " << run_action->GetVectorSize() << G4endl;

	for(int i=0; i<run_action->GetNumberOfRuns(); i++){
		G4cout << run_action->GetCrossVector(i) << " +/- " << run_action->GeterrVector(i) << " Barns || " << run_action->GetHitVector(i) << " +/- " 
		       << sqrt(run_action->GetHitVector(i)) << " hits || " << G4BestUnit(run_action->GetEnergyVector(i),"Energy") << G4endl;}
    
  delete runManager;
    
    

  return 0;
}

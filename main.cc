
#include "G4MTRunManager.hh"

#include "G4RunManager.hh"


#include "G4UImanager.hh"
#include "G4UserRunAction.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "QGSP_BIC_HP.hh"
#include "globals.hh"
#include "PhysicsList.hh"

int main(int argc, char** argv)
{
    G4UIExecutive* ui = 0;
  if ( argc == 1 ){
    ui = new G4UIExecutive(argc, argv);
  }

  G4Random::setTheEngine(new CLHEP::MixMaxRng);
  G4MTRunManager* runManager = new G4MTRunManager;

  // Set this true or false to test He3 or BF3, respectively:
  G4bool He3 = true;
  runManager->SetUserInitialization(new DetectorConstruction(He3));

  G4VModularPhysicsList* physicsList = new PhysicsList();
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  runManager->SetVerboseLevel(0);

  runManager->SetUserInitialization(new ActionInitialization(He3));

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( !ui ) {
    // batch mode - Apply macros directly 
    UImanager->ApplyCommand("/control/macroPath /home/codeyolson/he3_bf3_SourceFiles/macros");
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  } else {
    // Interactive Mode:
    //UImanager->ApplyCommand("/control/macroPath /reactorBay_sourceFiles");
    UImanager->ApplyCommand("/control/execute init.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;

}
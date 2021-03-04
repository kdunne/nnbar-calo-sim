//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "Analysis.hh"
#include "G4Types.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "PhysicsList.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include <G4ProductionCuts.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

//.........

std::vector<std::vector<G4double>> PMT_record;
std::vector<std::vector<G4double>> particle_gun_record;
std::vector<std::vector<G4double>> scint_record;
G4double event_number = 0;

std::ofstream PMT_outFile("./output/PMT_output_signal.txt");
std::ofstream Particle_outFile("./output/particle_output_signal.txt");
std::ofstream scint_outFile("./output/Scintillator_output_signal.txt");
std::ofstream SD_outFile("./output/All_SD_output_signal.txt");
std::ofstream TPC_outFile("./output/TPC_output_signal.txt");

std::ofstream Silicon_outFile("./output/Silicon_output_signal.txt");

//For the scintillator and lead glass position
std::ofstream Lead_glass_outFile("./output/Lead_glass_index.txt");
std::ofstream Scint_outFile("./output/Scint_index.txt");

/*** for gamma bkg 
std::ofstream PMT_outFile("./output/PMT_output_gamma_bkg.txt");
std::ofstream Particle_outFile("./output/particle_output_gamma_bkg.txt");
std::ofstream scint_outFile("./output/Scintillator_output_gamma_bkg.txt");
std::ofstream SD_outFile("./output/All_SD_output_gamma_bkg.txt");
***/ 

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " nnbar-calo-sim [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//.........

int main(int argc, char** argv)
{  

  scint_outFile << "Event_ID,group_ID,module_ID,Layer,Time,KE" << G4endl;
  Particle_outFile << "Event_ID,PID,Mass,KE,x,y,z,t,u,v,w"<< G4endl;
  SD_outFile<< "Event_ID,Track_ID,Parent_ID,Name,Proc,Volume,Group_ID,module_ID,index,x,y,z,t,KE,eDep,photons"<<G4endl;
  TPC_outFile << "Event_ID,module_ID,Layer,Track_ID,Parent_ID,Name,x,y,z,t,KE,eDep,electrons,trackl"<< G4endl;
  Lead_glass_outFile << "index,x,y,z" <<G4endl;
  Scint_outFile << "index,x,y,z" <<G4endl;
  Silicon_outFile << "Event_ID,layer,Track_ID,Parent_ID,Name,x,y,z,t,KE,eDep,trackl" <<G4endl;

  // Evaluate arguments
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Detect interactive mode (if no macro provided) and define UI session
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Optionally: choose a different Random engine...
  //
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the MT run manager
  //
#ifdef G4MULTITHREADED
  auto runManager = new G4MTRunManager;
  if ( nThreads > 0 ) { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#else
  auto runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);
 
  //G4VUserPhysicsList* physicsList = new PhysicsList();
  
 
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  opticalPhysics->SetWLSTimeProfile("delta");
  opticalPhysics->SetScintillationYieldFactor(1.0);
  opticalPhysics->SetScintillationExcitationRatio(0.0);
  opticalPhysics->SetMaxNumPhotonsPerStep(2000);
  opticalPhysics->SetMaxBetaChangePerStep(100.0);
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);
  physicsList->RegisterPhysics(opticalPhysics); 
   /***
  ***/
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(1.0*MeV, 10.0*TeV);

  

  runManager->SetUserInitialization(physicsList);

  auto actionInitialization = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);
  runManager->Initialize();
  
  // Initialize visualization
  auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

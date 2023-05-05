#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
//#include "Analysis.hh"
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
#include "QGSP_BERT.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

//.........

std::vector<std::vector<G4double>> PMT_record;
std::vector<std::vector<G4double>> particle_gun_record;
std::vector<std::vector<G4double>> scint_record;
G4double event_number = 0;
G4int run_number = 0;
// define sensitive volume readout :
std::ofstream Particle_outFile; 
std::ofstream Shield_outFile;

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " nnbar-veto-sandbox [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//.........

int main(int argc, char** argv)
{  

  // Evaluate arguments
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;

  #ifdef G4MULTITHREADED
    G4int nThreads = 4;
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
 
  G4VModularPhysicsList* physicsList = new PhysicsList();
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(30.0*eV, 10.0*TeV);
  runManager->SetUserInitialization(physicsList);
  
  auto actionInitialization = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  runManager->Initialize();
  
  // Initialize visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( macro.size() ) {
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

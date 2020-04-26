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

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager, ROOT specified in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms

  // Scintillator Bin Histos
  // 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  //
  //for(int i=0; i<10; i++) {
  //    G4String h1name = "Scint_Bin_" + i.c_str()
  analysisManager->CreateH1("Scint_Bin_1", "Num Scint Photons", 60000, 0, 60000); // 0
  analysisManager->CreateH1("Scint_Bin_2", "Num Scint Photons", 60000, 0, 60000); // 1
  analysisManager->CreateH1("Scint_Bin_3", "Num Scint Photons", 60000, 0, 60000); // 2
  analysisManager->CreateH1("Scint_Bin_4", "Num Scint Photons", 60000, 0, 60000); // 3
  analysisManager->CreateH1("Scint_Bin_5", "Num Scint Photons", 60000, 0, 60000); // 4
  analysisManager->CreateH1("Scint_Bin_6", "Num Scint Photons", 60000, 0, 60000); // 5
  analysisManager->CreateH1("Scint_Bin_7", "Num Scint Photons", 60000, 0, 60000); // 6
  analysisManager->CreateH1("Scint_Bin_8", "Num Scint Photons", 60000, 0, 60000); // 7
  analysisManager->CreateH1("Scint_Bin_9", "Num Scint Photons", 60000, 0, 60000); // 8
  analysisManager->CreateH1("Scint_Bin_10", "Num Scint Photons", 60000, 0, 60000); //9
  //    G4cout << "h1name: " << h1name << G4endl; 
  //}

  analysisManager->CreateH1("Eabs","Energy Deposited in Absorber", 500, 0., 500*MeV); //	10
  analysisManager->CreateH1("Egap","Energy Deposited in Scintillator", 100, 0., 1*MeV);	  // 	11
  analysisManager->CreateH1("Labs","Track Length in Absorber", 100, 0., 500.);  // 		12
  analysisManager->CreateH1("NumCerenkov", "Num Cerenkov Photons", 40000, 0, 40000); //      	13
  analysisManager->CreateH1("Lscint", "Track Length in Scintillator", 10, 0, 300.);        //   14
  analysisManager->CreateH1("TrackLength", "Total Track Length", 8000, 0, 800);  //          	15
  analysisManager->CreateH1("PhotonTime", "Cerenkov Photon Production", 50, 0, 3*ns); // 	16
  analysisManager->CreateH1("DecayTime", "Primary Decay Time", 50, 0, 3*ns); //         	17
  analysisManager->CreateH1("Range", "Primary Particle Range", 550, 0, 55); //  		18
  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax
  analysisManager->CreateH2("KinE","Kinetic Energy", 550, 0, 55, 350, 0, 350); // 		0 


  

}

//....

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//.....

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  
  //fileName = argv[1];
  // Open an output file
  //
  G4String fileName = "calo-sim";
  analysisManager->OpenFile(fileName);
}

//....

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....

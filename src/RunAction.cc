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
  analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV); // 0
  analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 800*MeV);	  // 1
  analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 500.);  // 2
  analysisManager->CreateH1("Cerenkov", "Num Cerenkov photons", 40000, 0, 40000); // 3
  analysisManager->CreateH1("Lscint", "trackL in scint", 10, 0, 300.);        // 4
  analysisManager->CreateH1("TrackLength", "Total Track Length", 8000, 0, 800);  // 5
  analysisManager->CreateH1("PhotonTime", "Cerenkov Photon Time Dist", 5000, 0, 5000*ns); // 6
  analysisManager->CreateH1("DecayTime", "Primary Decay Time", 5000, 0, 5000*ns); //7
  // create h2 energy dep time

  // Creating ntuple
  //
  analysisManager->CreateNtuple("calo", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Eabs");
  analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->CreateNtupleDColumn("Labs");
  analysisManager->CreateNtupleDColumn("Ceren");
  analysisManager->FinishNtuple();
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

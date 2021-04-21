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

#include <iostream>
#include <fstream>
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
//  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
//  analysisManager->SetNtupleMerging(true);

  // Book histograms

  G4cout << "Booking histograms " << G4endl;

  // 1-D Histos
  //analysisManager->CreateH1("NumCerenkov_2", "Num Cerenkov Photons Abs", 80, 0, 30000); 	//   	15
  analysisManager->CreateH1("NumChargedAbs"  ,"", 100, 0, 200 );   // 0
  analysisManager->CreateH1("NumChargedAct"  ,"", 100, 0, 100);   // 1
  analysisManager->CreateH1("NumPhoton"      ,"", 100, 0, 50000); // 2
  analysisManager->CreateH1("eDep"           ,"", 100, 0, 400);   // 3
  analysisManager->CreateH1("TotalNumCharged","", 100, 0,100 );   // 4
  analysisManager->CreateH1("FracNumCharged","", 100, 0, 1 );   //  5



  //analysisManager->CreateH1("NumCerenkov_2", "Num Cerenkov Photons Abs", 240, 0, 90000); 	//   	15
  // 2-D Histos
  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax
  //  analysisManager->CreateH2("A_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4


/***
  // eDep NTuple
  analysisManager->CreateNtuple("EnergyDeposit", "EnergyDeposit");  
  analysisManager->CreateNtupleDColumn("A0");	// 0
  analysisManager->CreateNtupleDColumn("A1");   // 1
  analysisManager->CreateNtupleDColumn("A2");   // 2
  analysisManager->CreateNtupleDColumn("A3");   // 3
  analysisManager->CreateNtupleDColumn("A4");   // 4
  analysisManager->CreateNtupleDColumn("A0");	// 5
  analysisManager->CreateNtupleDColumn("B1");   // 6
  analysisManager->CreateNtupleDColumn("B2");   // 7
  analysisManager->CreateNtupleDColumn("B3");   // 8
  analysisManager->CreateNtupleDColumn("B4");   // 9
  analysisManager->CreateNtupleDColumn("C0");	// 10
  analysisManager->CreateNtupleDColumn("C1");   // 11
  analysisManager->CreateNtupleDColumn("C2");   // 12
  analysisManager->CreateNtupleDColumn("C3");   // 13
  analysisManager->CreateNtupleDColumn("C4");   // 14
  analysisManager->CreateNtupleDColumn("D0");	// 15
  analysisManager->CreateNtupleDColumn("D1");   // 16
  analysisManager->CreateNtupleDColumn("D2");   // 17
  analysisManager->CreateNtupleDColumn("D3");   // 18
  analysisManager->CreateNtupleDColumn("D4");   // 19
  analysisManager->CreateNtupleDColumn("E0");	// 20
  analysisManager->CreateNtupleDColumn("E1");   // 21
  analysisManager->CreateNtupleDColumn("E2");   // 22
  analysisManager->CreateNtupleDColumn("E3");   // 23
  analysisManager->CreateNtupleDColumn("E4");   // 24
  analysisManager->FinishNtuple();
***/
}

//....

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//.....

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  G4String fileName = "../../output/lead-glass-sim.root";
//G4String fileName = "../../output/lead-glass-sim.csv";

  analysisManager->OpenFile(fileName);
}

//....

void RunAction::EndOfRunAction(const G4Run*)
{
  
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....

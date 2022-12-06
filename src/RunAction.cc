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
  analysisManager->SetVerboseLevel(1);

  // Book histograms
  G4cout << "Booking histograms " << G4endl;

  // 1-D Histos

  analysisManager->CreateH1("scintphotons_A", "photons", 2000, 0, 2000);       //      1
  analysisManager->CreateH1("scintphotons_B", "photons", 2000, 0, 2000);       //      2
  analysisManager->CreateH1("scintphotons_C", "photons", 2000, 0, 2000);       //      3
  analysisManager->CreateH1("scintphotons_D", "photons", 2000, 0, 2000);       //      4
  analysisManager->CreateH1("scintphotons_E", "photons", 2000, 0, 2000);       //      5
  analysisManager->CreateH1("scintphotons_F", "photons", 2000, 0, 2000);       //      6
  analysisManager->CreateH1("scintphotons_G", "photons", 2000, 0, 2000);       //      7
  analysisManager->CreateH1("scintphotons_H", "photons", 2000, 0, 2000);       //      8
  analysisManager->CreateH1("scintphotons_I", "photons", 2000, 0, 2000);       //      9
  analysisManager->CreateH1("scintphotons_J", "photons", 2000, 0, 2000);       //      10

  analysisManager->CreateH1("fiberAphotons_A", "photons", 700, 0, 2000);       //      11
  analysisManager->CreateH1("fiberAphotons_B", "photons", 700, 0, 2000);       //      12
  analysisManager->CreateH1("fiberAphotons_C", "photons", 700, 0, 2000);       //      13
  analysisManager->CreateH1("fiberAphotons_D", "photons", 700, 0, 2000);       //      14
  analysisManager->CreateH1("fiberAphotons_E", "photons", 700, 0, 2000);       //      15
  analysisManager->CreateH1("fiberAphotons_F", "photons", 700, 0, 2000);       //      16
  analysisManager->CreateH1("fiberAphotons_G", "photons", 700, 0, 2000);       //      17
  analysisManager->CreateH1("fiberAphotons_H", "photons", 700, 0, 2000);       //      18
  analysisManager->CreateH1("fiberAphotons_I", "photons", 700, 0, 2000);       //      19
  analysisManager->CreateH1("fiberAphotons_J", "photons", 700, 0, 2000);       //      20

  analysisManager->CreateH1("fiberBphotons_A", "photons", 700, 0, 2000);       //      21
  analysisManager->CreateH1("fiberBphotons_B", "photons", 700, 0, 2000);       //      22
  analysisManager->CreateH1("fiberBphotons_C", "photons", 700, 0, 2000);       //      23
  analysisManager->CreateH1("fiberBphotons_D", "photons", 700, 0, 2000);       //      24
  analysisManager->CreateH1("fiberBphotons_E", "photons", 700, 0, 2000);       //      25
  analysisManager->CreateH1("fiberBphotons_F", "photons", 700, 0, 2000);       //      26
  analysisManager->CreateH1("fiberBphotons_G", "photons", 700, 0, 2000);       //      27
  analysisManager->CreateH1("fiberBphotons_H", "photons", 700, 0, 2000);       //      28
  analysisManager->CreateH1("fiberBphotons_I", "photons", 700, 0, 2000);       //      29
  analysisManager->CreateH1("fiberBphotons_J", "photons", 700, 0, 2000);       //      30



/****
  analysisManager->CreateH1("Edep_A", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_B", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_C", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_D", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_E", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_F", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_G", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_H", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_I", "photons", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Edep_J", "photons", 700, 0, 700*MeV);       //      4


  // 2-D Histos
  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax

  analysisManager->CreateH2("A_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("B_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("C_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("D_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("E_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("F_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("G_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("H_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("I_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4
  analysisManager->CreateH2("J_eDepvscint", "eDepvscint", 60, 0, 60, 100, 0, 15000);			    // 4

  analysisManager->CreateH2("A_popFiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("B_popFiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("C_popFiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("D_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("E_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("F_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("G_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("H_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("I_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4
  analysisManager->CreateH2("J_popfiber", "fibervscint", 100, 0, 15000, 100, 0, 500);			    // 4


  // 3-D Histos
  // name, title, nxbins, xmin, xmas, nybins, ymin, ymax

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

// G4String fileName = "../../output/lead-glass-sim.root";
//G4String fileName = "../../output/lead-glass-sim.csv";
  G4String fileName = "../output/output.root";

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

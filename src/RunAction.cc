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
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // Book histograms

  G4cout << "Booking histograms " << G4endl;
// Scintillator Histograms Edeposited in individual sheets for low energy gamma backgrounds keV
/***
  analysisManager->CreateH1("Scint_Bin_1", "Energy Deposit in Scint Bin", 100, 0, 1000*keV); // 0
  analysisManager->CreateH1("Scint_Bin_2", "", 100, 0, 1000*keV); // 1
  analysisManager->CreateH1("Scint_Bin_3", "", 100, 0, 1000*keV); // 2
  analysisManager->CreateH1("Scint_Bin_4", "", 100, 0, 1000*keV); // 3
  analysisManager->CreateH1("Scint_Bin_5", "", 100, 0, 1000*keV); // 4
  analysisManager->CreateH1("Scint_Bin_6", "", 100, 0, 1000*keV); // 5
  analysisManager->CreateH1("Scint_Bin_7", "", 100, 0, 1000*keV); // 6
  analysisManager->CreateH1("Scint_Bin_8", "", 100, 0, 1000*keV); // 7
  analysisManager->CreateH1("Scint_Bin_9", "", 100, 0, 1000*keV); // 8
  analysisManager->CreateH1("Scint_Bin_10", "", 100, 0, 1000*keV); //9
***/


// Scintillator Histograms Edeposited in individual sheets for Primary Particles
//
  analysisManager->CreateH1("Scint_Bin_1", "Energy Deposit in Scint Bin", 100, 0, 100*MeV); // 0
  analysisManager->CreateH1("Scint_Bin_2", "", 100, 0, 100*MeV); // 1
  analysisManager->CreateH1("Scint_Bin_3", "", 100, 0, 100*MeV); // 2
  analysisManager->CreateH1("Scint_Bin_4", "", 100, 0, 100*MeV); // 3
  analysisManager->CreateH1("Scint_Bin_5", "", 100, 0, 100*MeV); // 4
  analysisManager->CreateH1("Scint_Bin_6", "", 100, 0, 100*MeV); // 5
  analysisManager->CreateH1("Scint_Bin_7", "", 100, 0, 100*MeV); // 6
  analysisManager->CreateH1("Scint_Bin_8", "", 100, 0, 100*MeV); // 7
  analysisManager->CreateH1("Scint_Bin_9", "", 100, 0, 100*MeV); // 8
  analysisManager->CreateH1("Scint_Bin_10", "", 100, 0, 100*MeV); //9

  // 1-D Histos
  analysisManager->CreateH1("NumCerenkov", "Num Cerenkov Photons", 80, 0, 10000); 		//   	10
  //analysisManager->CreateH1("NumCerenkov","NC", 80, 0, 500);
  analysisManager->CreateH1("PhotonTime", "Cerenkov Photon Production", 50, 0, 10*ns); 		// 	11
  analysisManager->CreateH1("DecayTime", "Primary Decay Time", 50, 0, 3*ns); 			//     	12 
  analysisManager->CreateH1("Range", "Primary Particle Range", 57, 0, 57); 			//	13
  //analysisManager->CreateH1("Range", "Primary Particle Range", 550, 0, 55); 			//	13
  analysisManager->CreateH1("EdepScint", "Energy Deposited in Scintillators", 50, 0, 500*MeV);  //      14
  analysisManager->CreateH1("EdepAbs", "Energy Deposited in Lead-glass", 50, 0, 500*MeV);       //      15
  analysisManager->CreateH1("EdepTube", "Energy Deposited in Vacuum Tube", 50, 0, 500*MeV);     //      16
  analysisManager->CreateH1("EdepElectron", "Energy Deposited by Electron", 50, 0, 10*MeV);      //      17
  analysisManager->CreateH1("RangeElectron", "Range of electron", 57, 0, 57);                   //      18

  // 2-D Histos
  // Two sets of histos
  // 500 MeV limit for signal
  // 4,10MeV limit for low energy background

  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax

  analysisManager->CreateH2("KinE","Kinetic Energy", 550, 0, 55, 350, 0, 350); 			             //	0 
  //analysisManager->CreateH2("KinE","Kinetic Energy", 57, 0, 57, 50, 0, 4*MeV); 		             // 0

  analysisManager->CreateH2("eDepvRange", "Energy Deposited", 57, 0, 57, 50, 0, 500*MeV );                   // 1
  //analysisManager->CreateH2("eDepvRange", "eDepvRange", 57, 0, 57, 50, 0, 10*MeV);                         // 1

  analysisManager->CreateH2("eDepvCerenkov", "Energy Deposited v Cerenkov", 80, 0, 10000, 50, 0, 500*MeV);   // 2
  //analysisManager->CreateH2("eDepvCerenkov", "eDepvCerenkov", 80, 0, 500, 50, 0, 3.5*MeV);                 // 2

  analysisManager->CreateH2("eDepvRangeElectron", "EDevRangeElectron", 57, 0, 57, 50, 0, 10*MeV);             // 3
  analysisManager->CreateH2("KinEelectron", "KinEelectron", 57, 0, 57, 50, 0, 10*MeV);                        // 4

  analysisManager->CreateH2("RangevCerenkov", "", 80, 0, 10000, 57, 0, 57); //5
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
  
   G4String fileName = "calo-sim.root";
  //G4String Dir = "../output/";
  //G4String fileName = Dir + "scint-calo-sim.root";
  //G4String fileName = Dir + "abs-calo-sim.root";
  //G4cout << "Printing to " << fileName << G4endl;

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

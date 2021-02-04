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
  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
//  analysisManager->SetNtupleMerging(true);

  // Book histograms

  G4cout << "Booking histograms " << G4endl;

  // 1-D Histos
  analysisManager->CreateH1("EdepAbs_A0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      1

/***
  analysisManager->CreateH1("EdepAbs_A1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_A9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      10


  analysisManager->CreateH1("EdepAbs_B0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      11
  analysisManager->CreateH1("EdepAbs_B1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_B9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      20

  analysisManager->CreateH1("EdepAbs_C0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      21
  analysisManager->CreateH1("EdepAbs_C1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_C9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      30


  analysisManager->CreateH1("EdepAbs_D0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      
  analysisManager->CreateH1("EdepAbs_D1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_D9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4

  analysisManager->CreateH1("EdepAbs_E0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_E9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4


  analysisManager->CreateH1("EdepAbs_F0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_F9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4

  analysisManager->CreateH1("EdepAbs_G0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_G9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4

  analysisManager->CreateH1("EdepAbs_H0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_H9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4

  analysisManager->CreateH1("EdepAbs_I0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_I9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4

  analysisManager->CreateH1("EdepAbs_J0", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J1", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J2", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J3", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J4", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J5", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J6", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J7", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J8", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
  analysisManager->CreateH1("EdepAbs_J9", "Energy Deposited in Lead-glass", 700, 0, 700*MeV);       //      4
***/


  analysisManager->CreateH1("NumCerenkov_Abs", "Num Cerenkov Photons Abs", 80, 0, 30000); 	//   	100

  analysisManager->CreateH1("PhotonTime", "Cerenkov Photon Production", 50, 0, 10*ns); 		// 	101
  analysisManager->CreateH1("DecayTime", "Primary Decay Time", 50, 0, 3*ns); 			//     	102 
  analysisManager->CreateH1("Range", "Primary Particle Range", 70, 0, 35); 			//	103
  analysisManager->CreateH1("EdepAbs", "Energy Deposited in Lead-glass", 50, 0, 700*MeV);       //      104
  analysisManager->CreateH1("Xpos", "Xpos", 100, -50, 50); 	                	        //      105	
  analysisManager->CreateH1("Ypos", "Ypos", 100, -50, 50); 	                		//      106	
  analysisManager->CreateH1("NumCerenkov_PMT", "Num Cerenkov Photons PMT", 80, 0, 60000); 	//   	107
  analysisManager->CreateH1("MoliereRadius", "MoliereRadius", 100, 0, 10); 			// 108
  analysisManager->CreateH1("EdepvRadius", "EdepvRadius", 50, 0, 50);				// 109
  analysisManager->CreateH1("EdepFrac", "EdepFrac", 100, 0, 1.1);				// 110
  analysisManager->CreateH1("gammaA_KE", "gammaA_KE", 130, 0, 650*MeV);                  // 111
  analysisManager->CreateH1("gammaB_KE", "gammaB_KE", 130, 0, 650*MeV) ;                 // 112
  analysisManager->CreateH1("angleDist", "angleDist", 72, 0, 360);			// 113
  analysisManager->CreateH1("endDist", "endDist", 60, 0, 6*m);				// 114

  //analysisManager->CreateH1("NumCerenkov_2", "Num Cerenkov Photons Abs", 80, 0, 30000); 	//   	115
  analysisManager->CreateH1("NumCerenkov_2", "Num Cerenkov Photons Abs", 240, 0, 60000); 	//   	116
  //analysisManager->CreateH1("NumCerenkov_2", "Num Cerenkov Photons Abs", 240, 0, 90000); 	//   	117

  // 2-D Histos
  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax
  analysisManager->CreateH2("RangevCerenkov", "Range v Cereknov", 35, 0, 35, 80, 0, 30000);      	    // 0
  analysisManager->CreateH2("eDepvCerenkov", "Energy Deposited v Cerenkov", 80, 0, 600*MeV, 80, 0, 60000);  // 1
  analysisManager->CreateH2("XY_Cerenkov", "XY Cerenkov", 100, -50, 50, 100, -50, 50);         	 	    // 2
  analysisManager->CreateH2("eDepvR", "eDepvR", 1000, 0, 10, 1000, 0, 10);				    // 3
  analysisManager->CreateH2("gammaKE", "gammaKE", 130, 0, 650*MeV, 130, 0, 650*MeV);			    // 4

  // 3-D Histos
  // name, title, nxbins, xmin, xmas, nybins, ymin, ymax
  analysisManager->CreateH3("lateralsize", "lateralsize", 100, -50, 50, 100, -50, 50, 100, 0, 1*keV); // 0

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

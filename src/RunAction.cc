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
  //analysisManager->SetNtupleMerging(true);

  // Book histograms

  G4cout << "Booking histograms " << G4endl;

  // 1-D Histos
  analysisManager->CreateH1("NumCerenkov_Abs", "Num Cerenkov Photons Abs", 80, 0, 30000); 	//   	0
  analysisManager->CreateH1("PhotonTime", "Cerenkov Photon Production", 50, 0, 10*ns); 		// 	1
  analysisManager->CreateH1("DecayTime", "Primary Decay Time", 50, 0, 3*ns); 			//     	2 
  analysisManager->CreateH1("Range", "Primary Particle Range", 35, 0, 35); 			//	3
  analysisManager->CreateH1("EdepAbs", "Energy Deposited in Lead-glass", 50, 0, 700*MeV);       //      4
  analysisManager->CreateH1("Xpos", "Xpos", 100, -50, 50); 	                	        //      5	
  analysisManager->CreateH1("Ypos", "Ypos", 100, -50, 50); 	                		//      6	
  analysisManager->CreateH1("NumCerenkov_PMT", "Num Cerenkov Photons PMT", 80, 0, 30000); 	//   	7
  analysisManager->CreateH1("MoliereRadius", "MoliereRadius", 100, 0, 10); 			// 8
  analysisManager->CreateH1("EdepvRadius", "EdepvRadius", 50, 0, 50);				// 9
  analysisManager->CreateH1("EdepFrac", "EdepFrac", 100, 0, 1);				// 10
  analysisManager->CreateH1("gammaA_KE", "gammaA_KE", 130, 0, 650*MeV);                  // 11
  analysisManager->CreateH1("gammaB_KE", "gammaB_KE", 130, 0, 650*MeV) ;                 // 12
  analysisManager->CreateH1("angleDist", "angleDist", 72, 0, 360);			// 13
  analysisManager->CreateH1("endDist", "endDist", 60, 0, 6*m);				// 14


  // 2-D Histos
  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax
  analysisManager->CreateH2("RangevCerenkov", "Range v Cereknov", 35, 0, 35, 80, 0, 30000);      	    // 0
  analysisManager->CreateH2("eDepvCerenkov", "Energy Deposited v Cerenkov", 80, 0, 30000, 50, 0, 250*MeV);  // 1
  analysisManager->CreateH2("XY_Cerenkov", "XY Cerenkov", 100, -50, 50, 100, -50, 50);         	 	    // 2
  analysisManager->CreateH2("eDepvR", "eDepvR", 1000, 0, 10, 1000, 0, 10);				    // 3
  analysisManager->CreateH2("gammaKE", "gammaKE", 130, 0, 650*MeV, 130, 0, 650*MeV);			    // 4

  // 3-D Histos
  // name, title, nxbins, xmin, xmas, nybins, ymin, ymax
  analysisManager->CreateH3("lateralsize", "lateralsize", 100, -50, 50, 100, -50, 50, 100, 0, 1*keV); // 0

  // eDep NTuple
  analysisManager->CreateNtuple("EnergyDeposit", "EnergyDeposit");  
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleDColumn("Radius");               //5
  analysisManager->CreateNtupleDColumn("eDep");               //7
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

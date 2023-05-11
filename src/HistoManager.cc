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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:fFactoryOn(false)
{
  	DefineCommands();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{ 
	// Creating a tree container.
	// This tree is associated to an output file.
	//
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	if ( ! fFactoryOn ) {
		//
		analysisManager->SetDefaultFileType("root");
		analysisManager->SetVerboseLevel(1);
		// Only merge in MT mode to avoid warning when running in Sequential mode
#ifdef G4MULTITHREADED
		analysisManager->SetNtupleMerging(true);
#endif
  		analysisManager->SetFileName(filename.c_str());
		G4bool fileOpen = analysisManager->OpenFile();
		if (! fileOpen) {
			G4cerr << "\n---> HistoManager::Book(): cannot open "
				<< analysisManager->GetFileName() << G4endl;
			return;
		}

		// create tree and branches
		analysisManager->CreateNtuple("nnbar", "NNbar simulation output");  
		analysisManager->CreateNtupleIColumn("p_evtno", p_evtno);
		analysisManager->CreateNtupleIColumn("p_pid", p_pid);
		analysisManager->CreateNtupleDColumn("p_m", p_m);
		analysisManager->CreateNtupleDColumn("p_q", p_q);
		analysisManager->CreateNtupleDColumn("p_ke", p_ke);
		analysisManager->CreateNtupleDColumn("p_x", p_x);
		analysisManager->CreateNtupleDColumn("p_y", p_y);
		analysisManager->CreateNtupleDColumn("p_z", p_z);
		analysisManager->CreateNtupleDColumn("p_t", p_t);
		analysisManager->CreateNtupleDColumn("p_px", p_px);
		analysisManager->CreateNtupleDColumn("p_py", p_py);
		analysisManager->CreateNtupleDColumn("p_pz", p_pz);
			
		analysisManager->CreateNtupleIColumn("det_pid", det_pid);
		analysisManager->CreateNtupleDColumn("det_t", det_t);
		analysisManager->CreateNtupleDColumn("det_ekin", det_ekin);
		analysisManager->CreateNtupleDColumn("det_x", det_x);
		analysisManager->CreateNtupleDColumn("det_y", det_y);
		analysisManager->CreateNtupleDColumn("det_z", det_z);
		analysisManager->CreateNtupleDColumn("det_px", det_px);
		analysisManager->CreateNtupleDColumn("det_py", det_py);
		analysisManager->CreateNtupleDColumn("det_pz", det_pz);
		analysisManager->FinishNtuple();  
		fFactoryOn = true;
	}

	G4cout << "\n----> Output file is open in "
		<< analysisManager->GetFileName() << "."
		<< analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{
	if (! fFactoryOn) return;

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

	G4cout << "\n----> Ntuple saved\n" << G4endl;
}

void HistoManager::ClearCryVectors()
{
	p_evtno.clear();
	p_pid.clear();
	p_m.clear();
	p_q.clear();
	p_ke.clear();
	p_x.clear();
	p_y.clear();
	p_z.clear();
	p_t.clear();
	p_px.clear();
	p_py.clear();
	p_pz.clear();
}


void HistoManager::ClearEventVectors()
{
	det_pid.clear();
	det_t.clear();
	det_ekin.clear();
	det_x.clear();
	det_y.clear();
	det_z.clear();
	det_px.clear();
	det_py.clear();
	det_pz.clear();
	
}

void HistoManager::FillCryVectors(G4int evtno, G4int pid, G4double m, G4double q,
			G4double ke, G4double x, G4double y, G4double z,
			G4double t, G4double px, G4double py, G4double pz)
{
	p_evtno.push_back(evtno);
	p_pid.push_back(pid);
	p_m.push_back(m);
	p_q.push_back(q);
	p_ke.push_back(ke);
	p_x.push_back(x);
	p_y.push_back(y);
	p_z.push_back(z);
	p_t.push_back(t);
	p_px.push_back(px);
	p_py.push_back(py);
	p_pz.push_back(pz);
}

void HistoManager::FillDetVectors(G4int pid, G4double t, G4double ekin, G4double x, G4double y, G4double z,
		G4double px, G4double py, G4double pz)
{
	det_pid.push_back(pid);
	det_t.push_back(t);
	det_ekin.push_back(ekin);
	det_x.push_back(x);
	det_y.push_back(y);
	det_z.push_back(z);
	det_px.push_back(px);
	det_py.push_back(py);
	det_pz.push_back(pz);
}

void HistoManager::FillTree()
{
	if (! fFactoryOn) return;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->AddNtupleRow(0);	
}

void HistoManager::DefineCommands()
{
	// Define /B5/generator command directory using generic messenger class
	fMessenger = new G4GenericMessenger(this, "/output/", "Output control");

	// randomizePrimary command
	auto& outputFileCmd = fMessenger->DeclareProperty("filename", filename);
	G4String guidance = "Path of output file.\n";
	outputFileCmd.SetGuidance(guidance);
	outputFileCmd.SetParameterName("outputFile", true);
	outputFileCmd.SetDefaultValue("Ana");
}


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
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:fFactoryOn(false)
{}

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
		G4bool fileOpen = analysisManager->OpenFile("ScintAna");
		if (! fileOpen) {
			G4cerr << "\n---> HistoManager::Book(): cannot open "
				<< analysisManager->GetFileName() << G4endl;
			return;
		}

		// create tree and branches
		analysisManager->CreateNtuple("nnbar", "NNbar simulation output");  
		analysisManager->CreateNtupleIColumn("p_evtno", p_evtno);
		analysisManager->CreateNtupleIColumn("p_parentid", p_parentid);
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
				
		analysisManager->CreateNtupleIColumn("scint_evtno", scint_evtno);
		analysisManager->CreateNtupleIColumn("scint_trackid", scint_trackid);
		analysisManager->CreateNtupleIColumn("scint_pid", scint_pid);
		analysisManager->CreateNtupleIColumn("scint_parentid", scint_parentid);
		analysisManager->CreateNtupleIColumn("scint_no", scint_no);
		analysisManager->CreateNtupleDColumn("scint_t", scint_t);
		analysisManager->CreateNtupleDColumn("scint_ke", scint_ke);
		analysisManager->CreateNtupleDColumn("scint_eDep", scint_eDep);
		analysisManager->CreateNtupleIColumn("scint_photons", scint_photons);
		analysisManager->CreateNtupleDColumn("scint_x", scint_x);
		analysisManager->CreateNtupleDColumn("scint_y", scint_y);
		analysisManager->CreateNtupleDColumn("scint_z", scint_z);
		analysisManager->CreateNtupleDColumn("scint_part_x", scint_part_x);
		analysisManager->CreateNtupleDColumn("scint_part_y", scint_part_y);
		analysisManager->CreateNtupleDColumn("scint_part_z", scint_part_z);
	
		analysisManager->CreateNtupleIColumn("fiber_evtno", fiber_evtno);
		analysisManager->CreateNtupleIColumn("fiber_trackid", fiber_trackid);
		analysisManager->CreateNtupleIColumn("fiber_parentid", fiber_parentid);
		analysisManager->CreateNtupleIColumn("fiber_procid", fiber_procid);
		analysisManager->CreateNtupleIColumn("fiber_no", fiber_no);
		analysisManager->CreateNtupleDColumn("fiber_t", fiber_t);
		analysisManager->CreateNtupleDColumn("fiber_x", fiber_x);
		analysisManager->CreateNtupleDColumn("fiber_y", fiber_y);
		analysisManager->CreateNtupleDColumn("fiber_z", fiber_z);
		analysisManager->CreateNtupleDColumn("fiber_part_x", fiber_part_x);
		analysisManager->CreateNtupleDColumn("fiber_part_y", fiber_part_y);
		analysisManager->CreateNtupleDColumn("fiber_part_z", fiber_part_z);
		
		analysisManager->CreateNtupleIColumn("sipm_evtno", sipm_evtno);
		analysisManager->CreateNtupleIColumn("sipm_trackid", sipm_trackid);
		analysisManager->CreateNtupleIColumn("sipm_parentid", sipm_parentid);
		analysisManager->CreateNtupleIColumn("sipm_no", sipm_no);
		analysisManager->CreateNtupleDColumn("sipm_t", sipm_t);
		analysisManager->CreateNtupleDColumn("sipm_ke", sipm_ke);
		analysisManager->CreateNtupleDColumn("sipm_x", sipm_x);
		analysisManager->CreateNtupleDColumn("sipm_y", sipm_y);
		analysisManager->CreateNtupleDColumn("sipm_z", sipm_z);

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

void HistoManager::ClearPVectors()
{
	p_evtno.clear();
	p_parentid.clear();
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
	scint_evtno.clear();
	scint_trackid.clear();
	scint_pid.clear();
	scint_parentid.clear();
	scint_no.clear();
	scint_t.clear();
	scint_ke.clear();
	scint_eDep.clear();
	scint_photons.clear();
	scint_x.clear();
	scint_y.clear();
	scint_z.clear();
	scint_part_x.clear();
	scint_part_y.clear();
	scint_part_z.clear();
	
	fiber_evtno.clear();
	fiber_trackid.clear();
	fiber_parentid.clear();
	fiber_procid.clear();
	fiber_no.clear();
	fiber_t.clear();
	fiber_x.clear();
	fiber_y.clear();
	fiber_z.clear();
	fiber_part_x.clear();
	fiber_part_y.clear();
	fiber_part_z.clear();
	
	sipm_evtno.clear();
	sipm_trackid.clear();
	sipm_parentid.clear();
	sipm_no.clear();
	sipm_t.clear();
	sipm_ke.clear();
	sipm_x.clear();
	sipm_y.clear();
	sipm_z.clear();

}


void HistoManager::FillPVectors(G4int evtno, G4int parentid, G4double m, G4double q,
			G4double ke, G4double x, G4double y, G4double z,
			G4double t, G4double px, G4double py, G4double pz)
{
	p_evtno.push_back(evtno);
	p_parentid.push_back(parentid);
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


void HistoManager::FillScintVectors(G4int evtno, G4int trackid, G4int pid,
		G4int parentid, G4int no, G4double t, G4double ke, G4double eDep,
		G4int photons, G4double x, G4double y, G4double z,
		G4double part_x, G4double part_y, G4double part_z)
{
	scint_evtno.push_back(evtno);
	scint_trackid.push_back(trackid);
	scint_pid.push_back(pid);
	scint_parentid.push_back(parentid);
	scint_no.push_back(no);
	scint_t.push_back(t);
	scint_ke.push_back(ke);
	scint_eDep.push_back(eDep);
	scint_photons.push_back(photons);
	scint_x.push_back(x);
	scint_y.push_back(y);
	scint_z.push_back(z);
	scint_part_x.push_back(part_x);
	scint_part_y.push_back(part_y);
	scint_part_z.push_back(part_z);

}


void HistoManager::FillFiberVectors(G4int evtno, G4int trackid, G4int parentid,
		G4int procid, G4int no, G4double t,
		G4double x, G4double y, G4double z,
		G4double part_x, G4double part_y, G4double part_z)
{
	fiber_evtno.push_back(evtno);
	fiber_trackid.push_back(trackid);
	fiber_parentid.push_back(parentid);
	fiber_procid.push_back(procid);
	fiber_no.push_back(no);
	fiber_t.push_back(t);
	fiber_x.push_back(x);
	fiber_y.push_back(y);
	fiber_z.push_back(z);
	fiber_part_x.push_back(part_x);
	fiber_part_y.push_back(part_y);
	fiber_part_z.push_back(part_z);

}


void HistoManager::FillSiPMVectors(G4int evtno, G4int trackid, G4int parentid,
		G4int no, G4double t, G4double ke,
		G4double x, G4double y, G4double z)
{
	sipm_evtno.push_back(evtno);
	sipm_trackid.push_back(trackid);
	sipm_parentid.push_back(parentid);
	sipm_no.push_back(no);
	sipm_t.push_back(t);
	sipm_ke.push_back(ke);
	sipm_x.push_back(x);
	sipm_y.push_back(y);
	sipm_z.push_back(z);
}


void HistoManager::FillTree()
{
	if (! fFactoryOn) return;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->AddNtupleRow(0);	
}


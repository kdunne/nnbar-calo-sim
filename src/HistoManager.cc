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

#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
	:fFactoryOn(false)
{
	DefineCommands();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
	delete G4AnalysisManager::Instance();
}

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
		analysisManager-> SetFirstNtupleId(0);
		// Only merge in MT mode to avoid warning when running in Sequential mode
#ifdef G4MULTITHREADED
		analysisManager->SetNtupleMerging(true);
#endif
  		analysisManager->SetFileName(filename.c_str());
		G4bool fileOpen = analysisManager->OpenFile();

		if (!fileOpen) {
			G4cerr << "\n---> HistoManager::Book(): cannot open "
				<< analysisManager->GetFileName() << G4endl;
			return;
		}

		analysisManager->Clear();
		this->ClearEventVectors();

		analysisManager->CreateNtuple("nnbar","NNBar main output");
			analysisManager->CreateNtupleIColumn("target_proc", target_proc);
		// analysisManager->FinishNtuple();

		auto scints = DetectorConstruction::GetPointer()->GetScintillatorNames();
		for(int sc = 0; sc < scints.size(); sc++)
		{
			/*
			// create tree and branches
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
			analysisManager->CreateNtupleDColumn("scint_ke", scint_ke);
			analysisManager->CreateNtupleDColumn("scint_eDep", scint_eDep);

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
			*/
				analysisManager->CreateNtupleIColumn(scints[sc]+"_scint_phot", scint_photons[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_scint_pdg", scint_pdg[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_scint_x", scint_x[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_scint_y", scint_y[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_scint_z", scint_z[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_scint_t", scint_t[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_holeA_phot", holeA_photons[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_holeB_phot", holeB_photons[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_fiberA_phot", fiberA_photons[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_fiberB_phot", fiberB_photons[sc]);

				analysisManager->CreateNtupleDColumn(scints[sc]+"_sipmA0_t", sipma0_t[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_sipmA1_t", sipma1_t[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_sipmB0_t", sipmb0_t[sc]);
				analysisManager->CreateNtupleDColumn(scints[sc]+"_sipmB1_t", sipmb1_t[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_sipmA0_phot", sipma0_phot[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_sipmA1_phot", sipma1_phot[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_sipmB0_phot", sipmb0_phot[sc]);
				analysisManager->CreateNtupleIColumn(scints[sc]+"_sipmB1_phot", sipmb1_phot[sc]);
		}
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

	// auto scints = DetectorConstruction::GetPointer()->GetScintillatorNames();
	// for(int sc = 0; sc < scints.size(); sc++)
	// {
	// 	G4cout << "HistoManager: scint_photons[" << scints[sc] << "] vector size: " << scint_photons[sc].size() << G4endl;
	// }

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

	fFactoryOn = false;

	G4cout << "\n----> Ntuple saved\n" << G4endl;
}

void HistoManager::ClearPVectors()
{
	/*
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
	*/
}


void HistoManager::ClearEventVectors()
{
	/*
	scint_evtno.clear();
	scint_trackid.clear();
	scint_pid.clear();
	scint_parentid.clear();
	scint_no.clear();
	scint_ke.clear();
	scint_eDep.clear();
	scintA_photons.clear();
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
	*/

	scint_photons.clear();
	holeA_photons.clear();
	holeB_photons.clear();
	fiberA_photons.clear();
	fiberB_photons.clear();

	sipma0_t.clear();
	sipma1_t.clear();
	sipmb0_t.clear();
	sipmb1_t.clear();
	sipma0_phot.clear();
	sipma1_phot.clear();
	sipmb0_phot.clear();
	sipmb1_phot.clear();

	scint_t.clear();
	scint_pdg.clear();
	scint_x.clear();
	scint_y.clear();
	scint_z.clear();

	// scint_z.clear();
	// scint_z.clear();

	target_proc.clear();

	for(int ii = 0; ii < DetectorConstruction::GetPointer()->GetScintillatorNames().size(); ii++)
	{
		std::vector<G4int> vecI = {};
		std::vector<G4double> vecD = {};

		scint_photons.push_back(vecI);
		scint_t.push_back(vecD);
		scint_pdg.push_back(vecI);
		scint_x.push_back(vecD);
		scint_y.push_back(vecD);
		scint_z.push_back(vecD);

		// scint_kEn.push_back(vecD);
		// scint_eDep.push_back(vecD);

		holeA_photons.push_back(vecI);
		holeB_photons.push_back(vecI);
		fiberA_photons.push_back(vecI);
		fiberB_photons.push_back(vecI);

		sipma0_t.push_back(vecD);
		sipma1_t.push_back(vecD);
		sipmb0_t.push_back(vecD);
		sipmb1_t.push_back(vecD);
		sipma0_phot.push_back(vecI);
		sipma1_phot.push_back(vecI);
		sipmb0_phot.push_back(vecI);
		sipmb1_phot.push_back(vecI);
	}
}


void HistoManager::FillPVectors(G4int evtno, G4int parentid, G4double m, G4double q,
			G4double ke, G4double x, G4double y, G4double z,
			G4double t, G4double px, G4double py, G4double pz)
{
	/*
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
	*/
}


void HistoManager::FillTargetVectors(G4int proc)
{
	target_proc.push_back(proc);
}

void HistoManager::FillScintVectors(G4int scint, G4int scintphotons,G4int holeAphotons, G4int holeBphotons, G4int fiberAphotons, G4int fiberBphotons, G4double hitTime, std::vector<G4double> hitPos, std::vector<G4int> hitPDG)
{
	auto scints = DetectorConstruction::GetPointer()->GetScintillatorNames();
	for(int sc = 0; sc < scints.size(); sc++)
	{
		if(sc!=scint) continue;

		if(scintphotons!=0)
			scint_photons[sc].push_back(scintphotons);
		if(holeAphotons!=0)
			holeA_photons[sc].push_back(holeAphotons);
		if(holeBphotons!=0)
			holeB_photons[sc].push_back(holeBphotons);
		if(fiberAphotons!=0)
			fiberA_photons[sc].push_back(fiberAphotons);
		if(fiberBphotons!=0)
			fiberB_photons[sc].push_back(fiberBphotons);

		scint_t[sc].push_back(hitTime);
		for(auto pdg : hitPDG)
			scint_pdg[sc].push_back(pdg);
		scint_x[sc].push_back(hitPos[0]);
		scint_y[sc].push_back(hitPos[1]);
		scint_z[sc].push_back(hitPos[2]);

		// for(auto eKin : hitkEn)
		// 	scint_kEn[sc].push_back(eKin);
		// for(auto eDep : hiteDep)
		// 	scint_eDep[sc].push_back(eDep);
	}
}

void HistoManager::FillFiberVectors(G4int evtno, G4int trackid, G4int parentid,
		G4int procid, G4int no, G4double t,
		G4double x, G4double y, G4double z,
		G4double part_x, G4double part_y, G4double part_z)
{
	/*
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
	*/
}


/*void HistoManager::FillSiPMVectors(const G4String& sipmname, G4double time, G4int phot)
{
	// sipm_evtno.push_back(evtno);
	// sipm_trackid.push_back(trackid);
	// sipm_parentid.push_back(parentid);
	// sipm_no.push_back(no);
	// sipm_t.push_back(t);
	// sipm_ke.push_back(ke);
	// sipm_x.push_back(x);
	// sipm_y.push_back(y);
	// sipm_z.push_back(z);

	// G4cout << "SIPM NAME:" << sipmname << " Time: " << time << " Phot: " << phot << G4endl;
	if(sipmname=="SiPMA0LV_A")
	{
		if (phot!=-1)
			sipma0_phot.push_back(phot);
		if (time!=0)
			sipma0_t.push_back(time);
	}
	else if (sipmname=="SiPMA1LV_A")
	{
		if (phot!=-1)
			sipma1_phot.push_back(phot);
		if (time!=0)
			sipma1_t.push_back(time);
	}
	else if (sipmname=="SiPMB0LV_A")
	{
		if (phot!=-1)
			sipmb0_phot.push_back(phot);
		if (time!=0)
			sipmb0_t.push_back(time);
	}
	else if (sipmname=="SiPMB1LV_A")
	{
		if (phot!=-1)
			sipmb1_phot.push_back(phot);
		if (time!=0)
			sipmb1_t.push_back(time);
	}
}*/

void HistoManager::FillSiPMVectors(const G4int scint, G4int photA0, G4int photA1, G4int photB0, G4int photB1)
{
	auto scints = DetectorConstruction::GetPointer()->GetScintillatorNames();
	for(int ii = 0; ii < scints.size(); ii++)
	{
		if(ii!=scint) continue;

		if(photA0!=0)
			sipma0_phot[ii].push_back(photA0);
		if(photA1!=0)
			sipma1_phot[ii].push_back(photA1);
		if(photB0!=0)
			sipmb0_phot[ii].push_back(photB0);
		if(photB1!=0)
			sipmb1_phot[ii].push_back(photB1);
	}
}

void HistoManager::FillSiPMVectors(const G4int scint, std::vector<G4double> tA0, std::vector<G4double> tA1, std::vector<G4double> tB0, std::vector<G4double> tB1)
{
	auto scints = DetectorConstruction::GetPointer()->GetScintillatorNames();
	for(int ii = 0; ii < scints.size(); ii++)
	{
		if(ii!=scint) continue;

		sipma0_t[ii] = tA0;
		sipma1_t[ii] = tA1;
		sipmb0_t[ii] = tB0;
		sipmb1_t[ii] = tB1;
	}
}

void HistoManager::FillTree() const
{
	if (! fFactoryOn) return;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	auto scints = DetectorConstruction::GetPointer()->GetScintillatorNames();

	G4int row = 0;

	analysisManager->AddNtupleRow(row);
	// row = 1;

	// for(int sc = 0; sc < scints.size(); sc++)
		// analysisManager->AddNtupleRow(sc+row);
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
	outputFileCmd.SetDefaultValue("ScintAna");
}

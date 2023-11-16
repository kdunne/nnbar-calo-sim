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
// #ifdef G4MULTITHREADED
// 		analysisManager->SetNtupleMerging(true);
// #endif
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

		analysisManager->CreateNtupleIColumn("sampling_trackid", sampling_trackid);
		analysisManager->CreateNtupleIColumn("sampling_pid", sampling_pid);
		analysisManager->CreateNtupleDColumn("sampling_t", sampling_t);
		analysisManager->CreateNtupleDColumn("sampling_ekin", sampling_ekin);
		analysisManager->CreateNtupleDColumn("sampling_x", sampling_x);
		analysisManager->CreateNtupleDColumn("sampling_y", sampling_y);
		analysisManager->CreateNtupleDColumn("sampling_z", sampling_z);
		analysisManager->CreateNtupleDColumn("sampling_px", sampling_px);
		analysisManager->CreateNtupleDColumn("sampling_py", sampling_py);
		analysisManager->CreateNtupleDColumn("sampling_pz", sampling_pz);
		
		analysisManager->CreateNtupleIColumn("cv_evtno", cv_evtno);
		analysisManager->CreateNtupleIColumn("cv_trackid", cv_trackid);
		analysisManager->CreateNtupleIColumn("cv_pid", cv_pid);
		analysisManager->CreateNtupleIColumn("cv_bar", cv_bar);
		analysisManager->CreateNtupleIColumn("cv_plane", cv_plane);
		analysisManager->CreateNtupleIColumn("cv_planedir", cv_planedir);
		//analysisManager->CreateNtupleSColumn("cv_name", cv_name);
		//analysisManager->CreateNtupleSColumn("cv_proc", cv_proc);
		//analysisManager->CreateNtupleSColumn("cv_vol_name", cv_vol_name);
		analysisManager->CreateNtupleDColumn("cv_t", cv_t);
		analysisManager->CreateNtupleDColumn("cv_ke", cv_ke);
		analysisManager->CreateNtupleDColumn("cv_eDep", cv_eDep);
		analysisManager->CreateNtupleDColumn("cv_tracklength", cv_tracklength);
		analysisManager->CreateNtupleDColumn("cv_x", cv_x);
		analysisManager->CreateNtupleDColumn("cv_y", cv_y);
		analysisManager->CreateNtupleDColumn("cv_z", cv_z);
		analysisManager->CreateNtupleDColumn("cv_px", cv_px);
		analysisManager->CreateNtupleDColumn("cv_py", cv_py);
		analysisManager->CreateNtupleDColumn("cv_pz", cv_pz);
		
		analysisManager->CreateNtupleIColumn("cvDigi_mul", cvDigi_mul);
		analysisManager->CreateNtupleIColumn("cvDigi_pid", cvDigi_pid);
		analysisManager->CreateNtupleIColumn("cvDigi_bar", cvDigi_bar);
		analysisManager->CreateNtupleIColumn("cvDigi_plane", cvDigi_plane);
		analysisManager->CreateNtupleDColumn("cvDigi_eDep", cvDigi_eDep);
		analysisManager->CreateNtupleDColumn("cvDigi_time", cvDigi_time);
		analysisManager->CreateNtupleDColumn("cvDigi_u", cvDigi_u);
		analysisManager->CreateNtupleDColumn("cvDigi_v", cvDigi_v);
		analysisManager->CreateNtupleDColumn("cvDigi_w", cvDigi_w);
		analysisManager->CreateNtupleDColumn("cvDigi_e1", cvDigi_e1);
		analysisManager->CreateNtupleDColumn("cvDigi_e2", cvDigi_e2);
		analysisManager->CreateNtupleDColumn("cvDigi_e3", cvDigi_e3);
		analysisManager->CreateNtupleDColumn("cvDigi_e4", cvDigi_e4);
		analysisManager->CreateNtupleDColumn("cvDigi_t1", cvDigi_t1);
		analysisManager->CreateNtupleDColumn("cvDigi_t2", cvDigi_t2);
		analysisManager->CreateNtupleDColumn("cvDigi_t3", cvDigi_t3);
		analysisManager->CreateNtupleDColumn("cvDigi_t4", cvDigi_t4);
		analysisManager->CreateNtupleDColumn("cvDigi_post", cvDigi_post);
		analysisManager->CreateNtupleDColumn("cvDigi_pose", cvDigi_pose);
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
	
	sampling_trackid.clear();
	sampling_pid.clear();
	sampling_t.clear();
	sampling_ekin.clear();
	sampling_x.clear();
	sampling_y.clear();
	sampling_z.clear();
	sampling_px.clear();
	sampling_py.clear();
	sampling_pz.clear();

	cv_evtno.clear();
	cv_trackid.clear();
	cv_pid.clear();
	cv_bar.clear();
	cv_plane.clear();
	cv_planedir.clear();
	//cv_name.clear();
	//cv_proc.clear();
	//cv_vol_name.clear();
	cv_t.clear();
	cv_ke.clear();
	cv_eDep.clear();
	cv_tracklength.clear();
	cv_x.clear();
	cv_y.clear();
	cv_z.clear();
	cv_px.clear();
	cv_py.clear();
	cv_pz.clear();
}

void HistoManager::ClearDigiVectors()
{
	cvDigi_mul.clear();
	cvDigi_mul.push_back(0);
	cvDigi_pid.clear();
	cvDigi_bar.clear();
	cvDigi_plane.clear();
	cvDigi_eDep.clear();
	cvDigi_time.clear();
	cvDigi_u.clear();
	cvDigi_v.clear();
	cvDigi_w.clear();
	cvDigi_e1.clear();
	cvDigi_e2.clear();
	cvDigi_e3.clear();
	cvDigi_e4.clear();
	cvDigi_t1.clear();
	cvDigi_t2.clear();
	cvDigi_t3.clear();
	cvDigi_t4.clear();
	cvDigi_post.clear();
	cvDigi_pose.clear();
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


void HistoManager::FillSamplingVectors(G4int trackid, G4int pid, G4double t, G4double ekin, G4double x, G4double y, G4double z,
		G4double px, G4double py, G4double pz)
{
	sampling_trackid.push_back(trackid);
	sampling_pid.push_back(pid);
	sampling_t.push_back(t);
	sampling_ekin.push_back(ekin);
	sampling_x.push_back(x);
	sampling_y.push_back(y);
	sampling_z.push_back(z);
	sampling_px.push_back(px);
	sampling_py.push_back(py);
	sampling_pz.push_back(pz);
}


void HistoManager::FillCVVectors(G4int evtno, G4int trackid, G4int pid, G4int bar, G4int plane,
		G4int planedir, G4double t, G4double ke, G4double eDep, 
		G4double tracklength, G4double x, G4double y, G4double z,
		G4double px, G4double py, G4double pz)
{
	cv_evtno.push_back(evtno);
	cv_trackid.push_back(trackid);
	cv_pid.push_back(pid);
	cv_bar.push_back(bar);
	cv_plane.push_back(plane);
	cv_planedir.push_back(planedir);
	cv_t.push_back(t);
	cv_ke.push_back(ke);
	cv_eDep.push_back(eDep);
	cv_tracklength.push_back(tracklength);
	cv_x.push_back(x);
	cv_y.push_back(y);
	cv_z.push_back(z);
	cv_px.push_back(px);
	cv_py.push_back(py);
	cv_pz.push_back(pz);

}


void HistoManager::FillCVDigiVectors(G4int pid, G4int bar, G4int plane, G4double eDep,
			G4double time, G4double u, G4double v, G4double w,
			G4double e1, G4double e2, G4double e3, G4double e4,
			G4double t1, G4double t2, G4double t3, G4double t4,
			G4double post, G4double pose)
{
	cvDigi_mul[0]++;
	cvDigi_pid.push_back(pid);
	cvDigi_bar.push_back(bar);
	cvDigi_plane.push_back(plane);
	cvDigi_eDep.push_back(eDep);
	cvDigi_time.push_back(time);
	cvDigi_u.push_back(u);
	cvDigi_v.push_back(v);
	cvDigi_w.push_back(w);
	cvDigi_e1.push_back(e1);
	cvDigi_e2.push_back(e2);
	cvDigi_e3.push_back(e3);
	cvDigi_e4.push_back(e4);
	cvDigi_t1.push_back(t1);
	cvDigi_t2.push_back(t2);
	cvDigi_t3.push_back(t3);
	cvDigi_t4.push_back(t4);
	cvDigi_post.push_back(post);
	cvDigi_pose.push_back(pose);

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


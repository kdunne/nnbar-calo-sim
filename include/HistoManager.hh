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
/// \file analysis/AnaEx01/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4AnalysisManager.hh"

class G4GenericMessenger;

class HistoManager
{
  public:
    HistoManager();
   ~HistoManager();
   
    void Book();
    void Save();

	void FillCryVectors(G4int evtno, G4int pid, G4double m, G4double q,
			G4double ke, G4double x, G4double y, G4double z,
			G4double t, G4double px, G4double py, G4double pz);
	void FillDetVectors(G4int pid, G4double t, G4double ekin, G4double x, G4double y, G4double z,
			G4double px, G4double py, G4double pz);
	void ClearCryVectors();
    void ClearEventVectors();

    void FillTree();
    
  private:
    void DefineCommands();
	
	G4String filename;
	G4GenericMessenger* fMessenger = nullptr;
	
	G4bool fFactoryOn;

	std::vector<G4int> p_evtno;
	std::vector<G4int> p_pid;
	std::vector<G4double> p_m;
	std::vector<G4double> p_q;
	std::vector<G4double> p_ke;
	std::vector<G4double> p_x;
	std::vector<G4double> p_y;
	std::vector<G4double> p_z;
	std::vector<G4double> p_t;
	std::vector<G4double> p_px;
	std::vector<G4double> p_py;
	std::vector<G4double> p_pz;
		
	std::vector<G4int> det_pid;
	std::vector<G4double> det_t;
	std::vector<G4double> det_ekin;
	std::vector<G4double> det_x;
	std::vector<G4double> det_y;
	std::vector<G4double> det_z;
   	std::vector<G4double> det_px;
	std::vector<G4double> det_py;
	std::vector<G4double> det_pz;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


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

#include "SteppingAction.hh"

#include <G4VSensitiveDetector.hh>

#include "EventAction.hh"
#include "Analysis.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

//....

G4double pmx, pmy, pmz, p_out, alfa_out, beta_out, gama_out, beal_out;

G4double a_cre, a_hb, a_tr, a_mf;
G4double tet, phi;

G4double  rad2deg = 180./3.141592654;



SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;

}

//....

SteppingAction::~SteppingAction()
{ 
}

//....

G4bool IsScintName(const G4String& name,const std::vector<G4String>& names)
{
	for(const auto& scintName : names)
		if(name=="ScintillatorLV_"+scintName)
			return true;

	return false;
}

G4bool IsLeftSDName(const G4String& name,const std::vector<G4String>& names)
{
	for(const auto& scintName : names)
		if(name=="SiPMA1LV_"+scintName || name=="SiPMB1LV_"+scintName)
			return true;

	return false;
}

G4bool IsRightSDName(const G4String& name,const std::vector<G4String>& names)
{
	for(const auto& scintName : names)
		if(name=="SiPMA0LV_"+scintName || name=="SiPMB0LV_"+scintName)
			return true;

	return false;
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{

	G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	if (eventNumber != fEventNumber) {
		fEventNumber = eventNumber;
		fScintillationCounter = 0;
		fCerenkovCounter = 0;
	}



	G4Track* track = step->GetTrack();
	G4int ID = track->GetTrackID();
	G4int ltime = track->GetLocalTime();

	const G4DynamicParticle* theParticle = track->GetDynamicParticle();

	G4ParticleDefinition *particleDef = track -> GetDefinition();
	G4String particleName =  particleDef -> GetParticleName();

	G4double eDep = step->GetTotalEnergyDeposit();
    G4int StepNumber  = step->GetTrack()->GetCurrentStepNumber();

	G4String material = track->GetMaterial()->GetName();

	G4int ParentID    = step->GetTrack()->GetParentID();
	G4double StepLen  = step->GetStepLength();
	G4double PartE    = step->GetTrack()->GetKineticEnergy();

	G4String ProcName;
	G4String CreProcName;

	const G4VProcess* Process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (Process) ProcName = Process->GetProcessName();

    const G4VProcess* CreProcess = track->GetCreatorProcess();
    if (CreProcess) CreProcName = CreProcess->GetProcessName();

	G4String proc;

	auto sdname = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();



//  if (particleName != "opticalphoton" && particleName != "e-")
// //if (sdname	== "TargetLV") {
// {
// G4cout  <<"EventNo: "<< eventNumber
//         <<" Particle: "<< particleName
//         <<" ID = "<<ID
//         <<" ParentID = "<<ParentID
//         <<" StepNumber = "<< StepNumber
//         <<" StepLen = " <<  G4BestUnit(StepLen, "Length")
//         <<" E = "<<  G4BestUnit(PartE, "Energy")
// //        PartE / MeV << " MeV"
//         <<" VolName = "<<sdname
//         <<" CreationProcName: " << CreProcName
//         <<" ProcName: " << ProcName
//         <<" IsFirstStepInVolume = "<< (step->IsFirstStepInVolume())
// //        <<" Step stst = "<< preStepStat
// //        <<" Track stat = "<< TrackStat
// //        <<" Tr_Coinc = "<< Tr_Coinc
// // 		  <<" inter = "<<inter
// <<G4endl;
// }



// 	if (ID != 1 && particleName!="proton") {
// 		proc = track->GetCreatorProcess()->GetProcessName();
//
// 		std::vector<G4String> sdnames = {"A"};
// //		auto sdname = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
// 		auto postGlobalTime = step->GetPostStepPoint()->GetGlobalTime();
//
// 		if (step->IsFirstStepInVolume())
// 		{
// 		}
//
// 		// G4cout << "ID: " << ID << " Event ID:" << eventNumber << " Particle name: " << particleName << " Material: " << material << " eDep: " << eDep << " Proccess name: " << proc << G4endl;
// 	} else {
// 		proc ="primary";
// 	}

//	if(particleName=="proton" && proc=="primary")

/*
	if(particleName=="proton" && proc=="primary")
	{
		auto det = step->GetPreStepPoint()->GetSensitiveDetector();

		if(det!=nullptr)*/



// if (
// 	(step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Target") &&
// 	(step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World") )
// {
//
// G4ThreeVector PartMom = theParticle->GetMomentum();
//
// pmx = PartMom.x();
// pmy = PartMom.y();
// pmz = PartMom.z();
//
// p_out = sqrt((pmx)*(pmx) + (pmy)*(pmy) + (pmz)*(pmz));
//
// alfa_out = pmx/p_out;
// beta_out = pmy/p_out;
// gama_out = pmz/p_out;
// beal_out = beta_out/alfa_out;
//
//
// tet=rad2deg*acos(gama_out);
// phi=rad2deg*atan(beal_out);
//
// //tet=rad2deg*asin(gama_out);
// //phi=rad2deg*asin(beta_out);
//
// //G4cout<<" tet = "<< tet<<" phi = "<<phi<<G4endl;
//
// if (abs(tet)<20. || abs(phi) > 1.44) {
// //	track->SetTrackStatus(fKillTrackAndSecondaries);
// track->SetTrackStatus(fStopAndKill);
// G4cout<<"  " << particleName<<" KILLED "<< G4endl;
// }
// // else G4cout<<"EventNo: "<< eventNumber
// //          <<" Particle: "<< particleName
// //          <<" ID = "<<ID
// //          <<" ProcName: " << ProcName
// //          <<" CreProcName = "<<CreProcName
// //          <<" tet = "<< tet
// //          <<" phi = "<<phi
// //          <<G4endl;
// }



} 


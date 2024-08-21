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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGenerator.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

namespace {G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(HistoManager *histo)
 : fHistoManager(histo)
{ 
	fPrimaryGenerator = new PrimaryGenerator();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{ 
	delete fPrimaryGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4AutoLock lock(&PrimaryGeneratorMutex);
	fPrimaryGenerator->GeneratePrimaryVertex(anEvent);
	
	// write particle properties to file
	fHistoManager->ClearPVectors();
	G4int evno = anEvent->GetEventID();
	G4int pid = fPrimaryGenerator->GetParticleDefinition1()->GetPDGEncoding();
	G4double mass = fPrimaryGenerator->GetParticleDefinition1()->GetPDGMass();
	G4double charge = fPrimaryGenerator->GetParticleDefinition1()->GetPDGCharge();
	G4double ke = fPrimaryGenerator->GetParticleEnergy1();
	G4double x = fPrimaryGenerator->GetParticlePosition().getX();
	G4double y = fPrimaryGenerator->GetParticlePosition().getY();
	G4double z = fPrimaryGenerator->GetParticlePosition().getZ();
	G4double t = 0.;
	G4double px = fPrimaryGenerator->GetParticleMomentumDirection1().getX();
	G4double py = fPrimaryGenerator->GetParticleMomentumDirection1().getY();
	G4double pz = fPrimaryGenerator->GetParticleMomentumDirection1().getZ();
	fHistoManager->FillPVectors(evno,pid,mass,charge,ke,x,y,z,t,px,py,pz);
		
	pid = fPrimaryGenerator->GetParticleDefinition2()->GetPDGEncoding();
	mass = fPrimaryGenerator->GetParticleDefinition2()->GetPDGMass();
	charge = fPrimaryGenerator->GetParticleDefinition2()->GetPDGCharge();
	ke = fPrimaryGenerator->GetParticleEnergy2();
	x = fPrimaryGenerator->GetParticlePosition().getX();
	y = fPrimaryGenerator->GetParticlePosition().getY();
	z = fPrimaryGenerator->GetParticlePosition().getZ();
	t = 0.;
	px = fPrimaryGenerator->GetParticleMomentumDirection2().getX();
	py = fPrimaryGenerator->GetParticleMomentumDirection2().getY();
	pz = fPrimaryGenerator->GetParticleMomentumDirection2().getZ();
	fHistoManager->FillPVectors(evno,pid,mass,charge,ke,x,y,z,t,px,py,pz);
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

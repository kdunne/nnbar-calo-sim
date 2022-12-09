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

#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

namespace {G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER;}
//.....

PrimaryGeneratorAction::PrimaryGeneratorAction(HistoManager *histo)
 : G4VUserPrimaryGeneratorAction(),fHistoManager(histo),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  // Hardcoded here for mu+ 50 MeV must be changed for different particle/momentum
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  fParticleGun->SetParticleDefinition(particleDefinition);
//  fParticleGun->SetParticlePosition(G4ThreeVector(15.,0.,0.));
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1.,0.,0.));
  fParticleGun->SetParticleEnergy(200*MeV);
}

//....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
	G4AutoLock lock(&PrimaryGeneratorMutex);
	// Set gun position
	fParticleGun->SetParticlePosition(G4ThreeVector(-20.*cm, 0.,0. ));

	fParticleGun->GeneratePrimaryVertex(anEvent);

	// write particle properties to file
	fHistoManager->ClearPVectors();
	G4int evno = anEvent->GetEventID();
	G4int pid = fParticleGun->GetParticleDefinition()->GetPDGEncoding();
	G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
	G4double charge = fParticleGun->GetParticleDefinition()->GetPDGCharge();
	G4double ke = fParticleGun->GetParticleEnergy();
	G4double x = fParticleGun->GetParticlePosition().getX();
	G4double y = fParticleGun->GetParticlePosition().getY();
	G4double z = fParticleGun->GetParticlePosition().getZ();
	G4double t = 0.;
	G4double px = fParticleGun->GetParticleMomentumDirection().getX();
	G4double py = fParticleGun->GetParticleMomentumDirection().getY();
	G4double pz = fParticleGun->GetParticleMomentumDirection().getZ();
	fHistoManager->FillPVectors(evno,pid,mass,charge,ke,x,y,z,t,px,py,pz);
	
}

//....


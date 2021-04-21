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

//.....

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  // Hardcoded here for mu+ 50 MeV must be changed for different particle/momentum
  auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
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


//  G4double worldZHalfLength = 35.*cm / 2.;
  G4double worldZHalfLength = 30.*cm / 2.;
  G4double gunPositionZ = -worldZHalfLength;
 
  //G4double worldZHalfLength = 42.*cm / 2.; // 35
  //G4double worldZHalfLength = 37.*cm / 2.; // 30
  //G4double worldZHalfLength = 32.*cm / 2.; // 25
  //G4double worldZHalfLength = 29.*cm/2;    // 22
//  G4double worldZHalfLength = 27.*cm/2.;   // 20
  //G4double worldZHalfLength = 25.*cm/2;      // 18
  //G4double worldZHalfLength = 22.*cm/2.;   // 15
  //G4double worldZHalfLength = 17.*cm/2.;   // 10
  //G4cout << "Gun at position " << (worldZHalfLength)/CLHEP::cm << " cm" <<G4endl;
  //G4double worldZHalfLength = 27.5*cm;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Set gun position
 // fParticleGun->SetParticlePosition(G4ThreeVector(-20.*cm, 0.,0. ));
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., gunPositionZ ));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
 

  std::cout << "Gun at Z Position: " << gunPositionZ/CLHEP::cm << std::endl;
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....


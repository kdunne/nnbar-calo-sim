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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#include "vector"
#include "RNGWrapper.hh"
#include "CRYSetup.h"
#include "PrimaryGeneratorMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"

using namespace std;


class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();

  public:
	void GeneratePrimaries(G4Event* anEvent);
	void InputCRY();
	void UpdateCRY(std::string* MessInput);
	void CRYFromFile(G4String newValue);
  
  // set methods
  void SetRandomFlag(G4bool value);

  private:
    G4ParticleGun*  fParticleGun; 

  private:
	G4GenericMessenger* fMessenger;
	std::vector<CRYParticle*>* vect; // vector of generated particles
	G4ParticleTable* particleTable;
	G4ParticleGun* particleGun;
	CRYGenerator* gen;
	G4int InputState;
	
	PrimaryGeneratorMessenger* gunMessenger;
};

//....

#endif

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
#include <iomanip>
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
using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction()
 :fParticleGun(nullptr),event_number(0)
{
  const char* inputfile = "setup.file";
  fParticleGun = new G4ParticleGun();

  std::ifstream inputFile;
  inputFile.open(inputfile, std::ios::in);
  
  char buffer[1000];

  if (inputFile.fail()) {
	  if (*inputfile != 0)  //....only complain if a filename was given
		  std::cerr << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << std::endl;
	  InputState = -1;
  }
  else {
	  std::string setupString("");
	  while (!inputFile.getline(buffer, 1000).eof()) {
		  setupString.append(buffer);
		  setupString.append(" ");
	  }

	CRYSetup* setup = new CRYSetup(setupString, "/mnt/e/cry_v1.7/data");

	gen = new CRYGenerator(setup);

	// set random number generator
	RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
	setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
	InputState = 0;

  }
  // create a vector to store the CRY particle properties
  vect = new std::vector<CRYParticle*>;
  particleTable = G4ParticleTable::GetParticleTable();
  gunMessenger = new PrimaryGeneratorMessenger(this);
  std::cerr << "Input state: " << InputState << std::endl;
  std::cerr << "particle table: " << particleTable << std::endl;
}

//....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
}

//....

void PrimaryGeneratorAction::InputCRY()
{
	InputState = 1;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::UpdateCRY(std::string* MessInput)
{
	CRYSetup* setup = new CRYSetup(*MessInput, "/mnt/e/cry_v1.7/data");

	gen = new CRYGenerator(setup);

	// set random number generator
	RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
	setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
	InputState = 0;

}

void PrimaryGeneratorAction::CRYFromFile(G4String newValue)
{
	// Read the cry input file
	std::ifstream inputFile;
	inputFile.open(newValue, std::ios::in);
	char buffer[1000];

	if (inputFile.fail()) {
		std::cerr << "Failed to open input file " << newValue << std::endl;
		std::cerr << "Make sure to define the cry library on the command line" << std::endl;
		InputState = -1;
	}
	else {
		
		std::string setupString("");
		while (!inputFile.getline(buffer, 1000).eof()) {
			setupString.append(buffer);
			setupString.append(" ");
		}

		CRYSetup* setup = new CRYSetup(setupString, "/mnt/e/cry_v1.7/data");

		gen = new CRYGenerator(setup);

		// set random number generator
		RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
		setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
		InputState = 0;
	}

	std::cerr << "Input state after CRYFromFile: " << InputState << std::endl;

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4String particleName;
	particleName = "mu-";
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu-");

	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleEnergy(500 * MeV);
	fParticleGun->SetParticlePosition(G4ThreeVector(0.* m, 0.* m, 2.0* m));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.*m, 0. * m, -1.0*m)); //(*vect)[j]->w())
	fParticleGun->SetParticleTime(1.*s);
	fParticleGun->GeneratePrimaryVertex(anEvent);


	event_ID.push_back(event_number);
	particle_ID.push_back(3); 
	particle_KE.push_back(50);
	particle_x.push_back(0.);
	particle_y.push_back(0.);
	particle_z.push_back(2.0);
	particle_momentum_x.push_back(0.);
	particle_momentum_y.push_back(0.);
	particle_momentum_z.push_back(1.0);
	particle_time.push_back(1.*s);
	
	event_number++;

}
	

//....


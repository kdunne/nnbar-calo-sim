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
// * technical work of the GEANEkin4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

namespace {G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER;}


PrimaryGenerator::PrimaryGenerator()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
	// Define setup geometry
	const G4double tgt_thickness = 2.5*mm;
	const G4double tgt_size = 1*cm;
	const G4double det_size = 5*cm;
	const G4double det_distance = 1*m;
	G4double phi_max = atan2(det_size/2.,det_distance);	// we will only generate particles in covered phi range.

	// Calculate reaction kinematics
	// Particle types
	partdef1 = G4ParticleTable::GetParticleTable()->FindParticle("proton");
	partdef2 = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
	
	// Particle mass and beam energy
	G4double m1 = partdef1->GetPDGMass();
	G4double m2 = partdef2->GetPDGMass();
	G4double m3 = m1;
	G4double m4 = m2; 
	G4double Ekin = 150.*MeV;

	// Incoming energy, momentum, etc
	G4double E1 = Ekin+m1;
	G4double E2 = m2;
	G4double p1 = sqrt(Ekin*Ekin+2*Ekin*m1);
	G4double p2 = 0;
	G4double beta = (p1+p2)/(E1+E2);
	G4double gamma = 1/(sqrt(1-beta*beta));
	G4double Ecm = sqrt(pow((E1+E2),2)-pow((p1+p2),2));
	G4double Ekincm = Ecm-m1-m2;

	// Outgoing energies etc in center-of-mass
	G4double Ekin3cm = (Ekincm/2)*(Ekincm+2*m4)/Ecm;
	G4double Ekin4cm = (Ekincm/2)*(Ekincm+2*m3)/Ecm;
	G4double E3cm = Ekin3cm+m3;
	G4double E4cm = Ekin4cm+m4;
	G4double p3cm = sqrt(Ekin3cm*Ekin3cm+2*Ekin3cm*m3);
	G4double p4cm = sqrt(Ekin4cm*Ekin4cm+2*Ekin4cm*m4);
	G4double beta3cm = p3cm/E3cm;
	G4double beta4cm = p4cm/E4cm;
	G4double pcm= sqrt(E3cm*E3cm-m3*m3);

	// Random ejectile angle in cm system
	G4double theta3cm=pi*G4UniformRand();
	// Ejectile
	G4double tantheta3 = sin(theta3cm)/(gamma*(cos(theta3cm)+beta/beta3cm));
	G4double theta3 = atan2(tantheta3,1);
	theta3 = (theta3<0) ? theta3+pi : theta3;
	G4double Ekin3 = (gamma-1)*m3+gamma*Ekin3cm+gamma*beta*pcm*cos(theta3cm);
	// Recoil
	G4double tantheta4 = sin(pi-theta3cm)/(gamma*(cos(pi-theta3cm)+beta/beta4cm));
	G4double theta4 = atan2(tantheta4,1);
	//theta4 = (theta4<0) ? theta4-pi : theta4;
	G4double Ekin4 = (gamma-1)*m4+gamma*Ekin4cm+gamma*beta*pcm*cos(pi-theta3cm);

	// Randomly generate phi within covered angular range
	G4double phi3 = 2*phi_max*G4UniformRand()-phi_max;
	G4double phi4 = phi3;
	G4double fiftyfifty = G4UniformRand();
	if(fiftyfifty<0.5){ phi3+=pi; }
	else{ phi4+=pi; }

	//particle 1 at vertex A
	G4ThreeVector mom1, mom2;
	mom1.setRThetaPhi(1,theta3,phi3);
	mom2.setRThetaPhi(1,theta4,phi4);
	
	// Create particles with energy and momentum from kinematics
	particle1 = new G4PrimaryParticle(partdef1);
	particle1->SetMomentumDirection(mom1);    
	particle1->SetKineticEnergy(Ekin3);
	
	particle2 = new G4PrimaryParticle(partdef2);
	particle2->SetMomentumDirection(mom2);    
	particle2->SetKineticEnergy(Ekin4);
		
	// Set random interaction point in target
	G4double x0=tgt_size*G4UniformRand()-tgt_size/2;
	G4double y0=tgt_size*G4UniformRand()-tgt_size/2;
	G4double z0=tgt_thickness*G4UniformRand()-tgt_thickness/2;
	
	position.set(x0,y0,z0-5*cm);
	G4double time = 0*s;
	
	G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);
	vertex->SetPrimary(particle1);
	vertex->SetPrimary(particle2);
	event->AddPrimaryVertex(vertex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

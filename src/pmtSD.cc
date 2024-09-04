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

#include "pmtSD.hh"
#include "NNbarHit.hh"

#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UserEventAction.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"


//......
pmtSD::pmtSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="pmtHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
    
}

//......
pmtSD::~pmtSD()
{}

//......
void pmtSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,
                                                             collectionName[0]);
}

//.....
G4bool pmtSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

    if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "PMT") return false;
    
    G4Track *theTrack = aStep->GetTrack();
    G4StepPoint* PostStep = aStep->GetPostStepPoint();

    //Get particle name
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    G4String particleName =  particleDef -> GetParticleName();
    
    // Get particle PDG code
    G4int pdg = particleDef->GetPDGEncoding();
    
    // Get unique track_id (in an event)
    G4int trackID = theTrack->GetTrackID();
    
    // Get Energy deposited
    G4double energyDeposit = aStep->GetTotalEnergyDeposit();
     
 
    // Get post-step position
    G4ThreeVector pos = PostStep->GetPosition();
//    G4double z = pos.getZ();
//    G4double x = pos.getX();
//    G4double y = pos.getY();

    // Get Vertex position
    G4ThreeVector vertex = theTrack->GetVertexPosition();
//    G4double vertZ = vertex.getZ();
//    G4double vertX = vertex.getX();
//    G4double vertY = vertex.getY();

    // Get Vertex Kinetic Energy
    G4double vertex_KE = theTrack->GetVertexKineticEnergy();

    // Read voxel indexes: i is the x index, k is the z index
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int k  = touchable->GetReplicaNumber(0);
    //G4int i  = touchable->GetReplicaNumber(2);
    //G4int j  = touchable->GetReplicaNumber(1);

    G4bool isLast = aStep->IsLastStepInVolume();

    // Get Global Time 
    G4double time = theTrack->GetGlobalTime() / CLHEP::ns;

    // Get Local Time
    G4double localTime = theTrack->GetLocalTime() / CLHEP::ns;


    // Should this be in Event Action?
    G4int photons = 0;
    G4ParticleDefinition* particle;
    if (particleName != "opticalphoton") {
        const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
        for (int j = 0; j < (*secondary).size(); j++) {
            particle = (*secondary)[j]->GetDefinition();
            if (particle->GetParticleName() == "opticalphoton" && (*secondary)[j]->GetCreatorProcess()->GetProcessName() == "Cerenkov") { photons++; } 
        }
    }
  
    
    // Get Process
    G4int parentID = 0;
    G4String proc = "";
    // Getting Process on parentID of primary causes seg fault
    if (trackID > 1){
	parentID = theTrack->GetParentID();
        proc = theTrack->GetCreatorProcess()->GetProcessName();
    } else {
        proc = "primary";
	parentID = 0;
    }

    if (proc=="Decay") {
        G4cout << "Killing particle " << particleName << G4endl;
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
	    
    // Get the pre-step kinetic energy
    G4double eKinPre = aStep -> GetPreStepPoint() -> GetKineticEnergy();
    // Get the post-step kinetic energy
    G4double eKinPost = aStep -> GetPostStepPoint() -> GetKineticEnergy();
    // Get the step average kinetic energy
    G4double eKinMean = (eKinPre + eKinPost) * 0.5;
        

/***
    if (trackID == 1 && isLast || eKinPost == 0) {
        const std::vector<const G4Track*>* lastSecondaries = aStep->GetSecondaryInCurrentStep();
        for (int j = 0; j < (*lastSecondaries).size(); j++) {
            particle = (*lastSecondaries)[j]->GetDefinition();
            G4String process = (*lastSecondaries)[j]->GetCreatorProcess()->GetProcessName();
            std::cout << "particle: " << particle->GetParticleName() << " process: " << process; 
        }

    } 
***/

    NNbarHit* detectorHit = new NNbarHit();

    // Particle Info
    detectorHit -> SetLocalTime(localTime);
    detectorHit -> SetParentID(parentID);
    detectorHit -> SetProcess(proc);
    detectorHit -> SetTime(time);
    detectorHit -> SetName(particleName);
    detectorHit -> SetTrackID(trackID);
    detectorHit -> SetIsLast(isLast);

    // Position Info
    detectorHit -> SetXID(k);
    detectorHit -> SetPos(pos);
    detectorHit -> SetVert(vertex);

    // Energy Info
    detectorHit -> SetEDep(energyDeposit);
    detectorHit -> SetVertexKE(vertex_KE);
    detectorHit -> SetKinEn(eKinPost);

    HitsCollection -> insert(detectorHit);

//	G4cout << "Replica: "       << k << G4endl;
//	G4cout << "tracklength: "   << tracklength << G4endl;
//	G4cout << "energyDeposit: " << energyDeposit << G4endl;
//	G4cout << "eKinPost: "      << eKinPost << G4endl;

    return true;
}

//......
void pmtSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}


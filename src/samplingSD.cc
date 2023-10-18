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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "samplingSD.hh"
//#include "TrackInformation.hh"

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
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

//.....
samplingSD::samplingSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname="samplingHitCollection";
    collectionName.insert(HCname);
    HitsCollection = NULL;
    sensitiveDetectorName = name;
}

//.....
samplingSD::~samplingSD()
{}

//.....
void samplingSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);
}

//.....
G4bool samplingSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

    //std::cout<< " shiled hit class SD " << std::endl;

    //if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "plane") return false;
    
    // Get Direction
    G4Track * theTrack = aStep  ->  GetTrack();
 
    // Get step length  
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    G4StepPoint* PostStep = aStep->GetPostStepPoint();

    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    G4String namePre = volumePre->GetName();
 
    G4TouchableHandle touchPostStep = PostStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePost = touchPostStep->GetVolume();
    G4String namePost = volumePost->GetName();
   
  	//TrackInformation* trackInfo = (TrackInformation*)(theTrack->GetUserInformation());
  	//trackInfo->SetTrackingStatus(1);
    G4int trackID =  theTrack->GetTrackID();
    G4int parentID =  theTrack->GetParentID();
	//if(){
	if(PostStep->GetStepStatus()==fGeomBoundary){
	
			//G4cout << trackID << G4endl;
			//G4cout << theTrack->GetCurrentStepNumber() << "," << aStep->IsFirstStepInVolume() << "," << aStep->IsLastStepInVolume() << "," << namePre << "," << namePost << G4endl;
		//}

		G4int pid = theTrack->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
		// Position
		G4ThreeVector pos1 = PreStep->GetPosition();
		G4ThreeVector pos2 = PostStep->GetPosition();
		G4double x = ((pos1+pos2)/2.).getX();
		G4double y = ((pos1+pos2)/2.).getY();
		G4double z = ((pos1+pos2)/2.).getZ();

		//G4cout << x << "," << y << "," << z << G4endl;
		G4ThreeVector momentum = PreStep->GetMomentumDirection();
		G4double px = momentum.getX();
		G4double py = momentum.getY();
		G4double pz = momentum.getZ();
		//G4cout << px << "," << py << "," << pz << G4endl;

		// Get Time
		G4double time = theTrack->GetGlobalTime() / CLHEP::ns;
		//G4cout << time << G4endl;
		
		// Get the kinetic energy
		G4double ekin = (PreStep->GetKineticEnergy()+PostStep->GetKineticEnergy())/2.;
  
		G4String name = theTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();

		NNbarHit* samplingHit = new NNbarHit();

		// Make this kinetic energy and position
		samplingHit -> SetTrackID(trackID);
		samplingHit -> SetPID(pid);
		samplingHit -> SetTime(time);
		samplingHit -> SetName(name);
		samplingHit -> SetKinEn(ekin);
		samplingHit -> SetPosX(x);
		samplingHit -> SetPosY(y);
		samplingHit -> SetPosZ(z);
		samplingHit -> SetPX(px);
		samplingHit -> SetPY(py);
		samplingHit -> SetPZ(pz);
		HitsCollection -> insert(samplingHit);
	}
    //}
    return true;
}

//......
void samplingSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}


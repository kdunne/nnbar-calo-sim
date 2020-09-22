#include "ScintillatorSD.hh"

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


//.....
ScintillatorSD::ScintillatorSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="ScintillatorHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
    
}

//.....
ScintillatorSD::~ScintillatorSD()
{}

//.....
void ScintillatorSD::Initialize(G4HCofThisEvent*) 
{HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);}

//.....
G4bool ScintillatorSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
    if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "Layer") return false;

    // Get Direction
    G4Track * theTrack = aStep  ->  GetTrack();
    G4ThreeVector stepDelta = aStep->GetDeltaPosition();
    G4double direction = stepDelta.getZ();

    //Get particle name
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    G4String particleName =  particleDef -> GetParticleName();
    
    // Get particle PDG code
    G4int pdg = particleDef ->GetPDGEncoding();
    
    // Get unique track_id (in an event)
    G4int trackID = theTrack -> GetTrackID();
   
    // Get Energy deposited
    //G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
    G4double energyDeposit = aStep->GetPreStepPoint()->GetKineticEnergy() - aStep->GetPostStepPoint()->GetKineticEnergy();
 
    // Get step length  
    G4double DX = aStep -> GetStepLength();
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    G4StepPoint* PostStep = aStep->GetPostStepPoint();    

    // Position z=0 at center of detector
    G4ThreeVector pos = PostStep->GetPosition();
    //G4ThreeVector pos = PreStep->GetPosition();
    G4double z = pos.getZ();

    G4ThreeVector vertex = theTrack->GetVertexPosition();
    G4double origin = vertex.getZ();
    G4double tracklength = z - origin;

    // Read voxel indexes: i is the x index, k is the z index
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int k  = touchable->GetReplicaNumber(0);
	G4int origin_replica = theTrack->GetOriginTouchable()->GetReplicaNumber(0); // ** I added this here !!!
    //G4int i  = touchable->GetReplicaNumber(2);
    //G4int j  = touchable->GetReplicaNumber(1);
  
    // Get Time
    G4double time = theTrack->GetGlobalTime() / CLHEP::ns;

    // Get Local Time
    G4double localTime = theTrack->GetLocalTime() / CLHEP::ns;

    // Get Name

    G4String name = theTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    G4String namePre = volumePre->GetName();

    G4int parentID = 0;
    G4String proc = ""; 

    // Get Process
    if (trackID > 1){
        parentID = theTrack->GetParentID();
        //proc = theTrack->GetCreatorProcess()->GetProcessName(); // disabled because of the error! 
    } else {
        proc = "primary";
	    parentID = 0;
        if (aStep->GetPreStepPoint()->GetKineticEnergy() == 0){theTrack->SetTrackStatus(fKillTrackAndSecondaries);}
    }

    if (proc=="Decay"){
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        return false; 
    }

	// Get the pre-step kinetic energy
    G4double eKinPre = aStep -> GetPreStepPoint() -> GetKineticEnergy();
    // Get the post-step kinetic energy
    G4double eKinPost = aStep -> GetPostStepPoint() -> GetKineticEnergy();

	NNbarHit* detectorHit = new NNbarHit();

    // Make this kinetic energy and position
    detectorHit -> SetLocalTime(localTime);
    detectorHit -> SetParentID(parentID);
    detectorHit -> SetProcess(proc);
    detectorHit -> SetTime(time);
    detectorHit -> SetName(name);

    detectorHit -> SetTrackID(trackID);
    detectorHit -> SetXID(k);
    detectorHit -> SetPosZ(tracklength);
    detectorHit -> SetEDep(energyDeposit);
    detectorHit -> SetKinEn(eKinPost);

	detectorHit->SetOrigin(origin_replica);


    HitsCollection -> insert(detectorHit);

    return true;
}

//......
void ScintillatorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	//std::cout << "** 6.5 Scintillator Event END --- " << std::endl;

    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}


#include "GenericSD.hh"
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
GenericSD::GenericSD(G4String name):
	G4VSensitiveDetector(name)
{
	collectionName.insert(name+"_Collection");
	HitsCollection = NULL;
	sensitiveDetectorName = name;
}

//......
GenericSD::~GenericSD()
{}

//......
void GenericSD::Initialize(G4HCofThisEvent*)
{
	HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);
}

//.....
G4bool GenericSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

	G4Track *theTrack = aStep->GetTrack();
	G4StepPoint* PreStep = aStep->GetPreStepPoint();

	//Get particle name
	G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
	G4String particleName =  particleDef -> GetParticleName();
	//const G4DynamicParticle *dynParticle = theTrack->GetDynamicParticle();

	G4int pdg = particleDef->GetPDGEncoding(); // Get particle PDG code
	G4int trackID = theTrack->GetTrackID(); // Get unique track_id (in an event)
	G4double eKin = theTrack->GetKineticEnergy(); // Get kinetic energy
	G4double eDep = aStep->GetTotalEnergyDeposit(); // Get deposited energy
	G4ThreeVector pos = PreStep->GetPosition(); // Get pre-step position
	G4ThreeVector vertex = theTrack->GetVertexPosition(); // Get Vertex position
	G4double vertex_KE = theTrack->GetVertexKineticEnergy(); // Get Vertex Kinetic Energy

	// Read voxel indexes: i is the x index, k is the z index
	const G4VTouchable* touchable = PreStep->GetTouchable();
	G4int k  = touchable->GetReplicaNumber(0);
	G4String detname  = touchable->GetVolume()->GetName();

	G4bool isLast = aStep->IsLastStepInVolume();
	G4double time = theTrack->GetGlobalTime() / CLHEP::ns; // Get Global Time 
	G4double localTime = theTrack->GetLocalTime() / CLHEP::ns; // Get Local Time

	// Should this be in Event Action?
	G4int photons = 0;
	G4ParticleDefinition* particle;
	if (particleName != "opticalphoton") {
		const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
		for (int j = 0; j < (*secondary).size(); j++) {
			particle = (*secondary)[j]->GetDefinition();
			if (particle->GetParticleName() == "opticalphoton" && (*secondary)[j]->GetCreatorProcess()->GetProcessName() == "Scintillation") { photons++; } 
		}
	}

	// Get Process
	G4int parentID = 0;
	G4String proc = "";
	// Getting Process on parentID of primary causes seg fault
	if(trackID > 1){
		parentID = theTrack->GetParentID();
		proc = theTrack->GetCreatorProcess()->GetProcessName();
	}
	else {
		proc = "primary";
		parentID = 0;
	}

	if(proc=="Decay") {
		G4cout << "Killing particle " << particleName << G4endl;
		theTrack->SetTrackStatus(fKillTrackAndSecondaries);
	}
	
	NNbarHit* detectorHit = new NNbarHit();

	// Particle Info
	detectorHit -> SetLocalTime(localTime);
	detectorHit -> SetParentID(parentID);
	detectorHit -> SetProcess(proc);
	detectorHit -> SetTime(time);
	detectorHit -> SetName(particleName);
	detectorHit -> SetParticleID(pdg);
	detectorHit -> SetTrackID(trackID);
	detectorHit -> SetIsLast(isLast);

	// Position Info
	detectorHit -> SetXID(k);
	detectorHit -> SetPos(pos);
	detectorHit -> SetVert(vertex);
	detectorHit -> SetDetName(detname);

	// Energy Info
	detectorHit -> SetEDep(eDep);
	detectorHit -> SetVertexKE(vertex_KE);
	detectorHit -> SetKinEn(eKin);
	detectorHit -> SetPhotons(photons);

	HitsCollection -> insert(detectorHit);

//	G4cout << "Processing event in " <<touchable->GetVolume()->GetName() << G4endl;
//	G4cout << "Detector " << HitsCollection->GetSDname() << " Collection " << HitsCollection->GetName() << G4endl;
//	G4cout << "Collection size " << HitsCollection->entries() << G4endl;
	return true;
}

void GenericSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	static G4int HCID = -1;
	if(HCID < 0)
	{
		HCID = GetCollectionID(0);
	}
	HCE -> AddHitsCollection(HCID,HitsCollection);
}

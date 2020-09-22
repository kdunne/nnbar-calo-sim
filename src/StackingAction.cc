
#include "StackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4LogicalVolume.hh"
#include "RunAction.hh"
#include <G4VPhysicalVolume.hh>
#include <vector>

StackingAction::StackingAction()
	: cerenkovCounter_all(0), cerenkovCounter_pri(0), scintCounter(0)
{}

StackingAction::~StackingAction()
{}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
	if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
	{ // particle is optical photon
		if (aTrack->GetParentID()>0)
		{ // particle is secondary
			if (aTrack->GetCreatorProcess()->G4VProcess::GetProcessName()=="Cerenkov" && 
				aTrack->GetOriginTouchable()->GetVolume(0)->G4VPhysicalVolume::GetName() == "AbsoPV"){
				cerenkovCounter_all ++;  
				if (aTrack->GetParentID() == 1) { cerenkovCounter_pri++; }
			}

			/***
			if (aTrack->GetCreatorProcess()->G4VProcess::GetProcessName() == "Scintillation" 
				&& aTrack->GetOriginTouchable()->GetVolume(0)->G4VPhysicalVolume::GetName() == "Layer") {
				auto index = aTrack->GetTouchable()->GetReplicaNumber(0);

				//G4cout << "#### can I get origin layer " << aTrack->GetOriginTouchable()->GetReplicaNumber(0) << G4endl;
				scint_layer_photon[index] = scint_layer_photon[index]+1;
			} ***/

		}
	}
	return fUrgent;
}


void StackingAction::NewStage()
{
	std::cout << "Number of Cerenkov photons: " << cerenkovCounter_pri << " all : " << cerenkovCounter_all<< std::endl;
	lead_glass_photon_all.push_back(cerenkovCounter_all);
	lead_glass_photon_pri.push_back(cerenkovCounter_pri);
	//scint_photon_collection.push_back(scint_layer_photon);
}

void StackingAction::PrepareNewEvent()
{
	cerenkovCounter_pri = 0; cerenkovCounter_all = 0;
	scintCounter = 0; // this is no longer useful 
	std::fill(scint_layer_photon.begin(), scint_layer_photon.end(), 0); // reset the collection array
}

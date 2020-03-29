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

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "Analysis.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

//....

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;

  //auto analysisManager = G4AnalysisManager::Instance();
  //auto hist = analysisManager->CreateH1("TimeDist", "Time Dist", 5000, 0, 5000);
  // analysisManager->OpenFile("timedist.csv");

}

//....

SteppingAction::~SteppingAction()
{ 
  //auto analysisManager = G4AnalysisManager::Instance(); 
  //analysisManager->Write();
  //analysisManager->CloseFile();
}

//....

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventNumber = G4RunManager::GetRunManager()->
                                              GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
     fEventNumber = eventNumber;
     fScintillationCounter = 0;
     fCerenkovCounter = 0;
  }

  G4Track* track = step->GetTrack();

  G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();
  if (track->GetTrackID() > 1) {
    if (ParticleName != "e-" && ParticleName != "opticalphoton") {
      G4cout << "Killing particle: " << ParticleName << G4endl;
      track->SetTrackStatus(fStopAndKill);
    }
  }


  if (ParticleName == "opticalphoton" && track->GetCreatorProcess()->GetProcessName()=="Cerenkov") {
  
    
 // G4double time = track->GetGlobalTime();
 // G4double loctime = track->GetLocalTime();

 // G4cout << "Local TIme :" << loctime << G4endl;
 // G4cout << "Globat Time :" << time << G4endl;

 // auto analysisManager = G4AnalysisManager::Instance();
 // analysisManager->FillH1(6, time);
  }

//  return;
  

  const std::vector<const G4Track*>* secondaries =
                                            step->GetSecondaryInCurrentStep();

  if (secondaries->size()>0) {
     for(unsigned int i=0; i<secondaries->size(); ++i) {
        if (secondaries->at(i)->GetParentID()>0) {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4OpticalPhoton::OpticalPhotonDefinition()){
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()== "Scintillation"){		
		fScintillationCounter++;
  		G4StepPoint* preStepPoint = step->GetPreStepPoint();
  		G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
 // 		G4int copyNo = theTouchable->GetReplicaNumber(); 
//		G4cout << "copyNumber : " << copyNo << G4endl;
	      }


	      if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Cerenkov")fCerenkovCounter++;
           }
        }
     }
  }

  else{
    G4int TrackID = track->GetTrackID(); 
    if (TrackID == 1 && track->GetKineticEnergy() == 0) {
      //G4cout << "stopping position" << track->GetPosition()  << G4endl;
    }
  }
}

//....

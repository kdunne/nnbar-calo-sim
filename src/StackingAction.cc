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

#include "StackingAction.hh"
#include "Analysis.hh"
#include "G4RunManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

//....

StackingAction::StackingAction()
  : G4UserStackingAction(),
    fScintillationCounter(0), fCerenkovCounter(0), felectroncounter(0), fposcounter(0), fgammacounter(0), fnumucounter(0), fantinuecounter(0),
    felectronmom(0), fgammamom(0), fnumumom(0), fantinuemom(0)
{}

//....

StackingAction::~StackingAction()
{}

//....

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{



  // muon decay products specified here for sanity checking the results

  G4ParticleDefinition* pos = G4Positron::PositronDefinition();
  G4ParticleDefinition* y = G4Gamma::GammaDefinition();
  G4ParticleDefinition* antinu_mu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4ParticleDefinition* nu_e = G4NeutrinoE::NeutrinoEDefinition(); 
  G4ParticleDefinition* opticalphoton = G4OpticalPhoton::OpticalPhotonDefinition();	
  G4ParticleDefinition* e = G4Electron::ElectronDefinition();

  //G4String proc = aTrack->GetCreatorProcess()->G4VProcess::GetProcessName(); 
  G4ParticleDefinition* pid = aTrack->GetDefinition();

  if(pid == opticalphoton)
  { // particle is optical photon
          if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
              fScintillationCounter++;
          if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov") {
          

	      // Time Create 
              G4double time = aTrack->GetGlobalTime();
              //G4cout << " Time ------ " << time << G4endl;
	      auto analysisManager = G4AnalysisManager::Instance();
              analysisManager->FillH1(6, time);
         

              //G4double time = aTrack->GetGlobalTime();
              //G4double loctime = aTrack->GetLocalTime();

              //G4cout << "Local Time : " << loctime/CLHEP::ns << G4endl;
              G4cout << "Global Time : " << time/CLHEP::ns << G4endl;
	      if(aTrack->GetParentID() > 1 && time/CLHEP::ns >=3 ) {
	     
	      G4cout << "Cerenkov Parent Track ID : " << aTrack->GetParentID() << G4endl; 
	      }//G4cout << "Volume : " << aTrack->GetVolume() << G4endl;

	      // Parent Particle
              //G4Track* pParent = G4Track(aTrack->GetParentID());         
	      //analysisManager->FillH1(8, position);


              fCerenkovCounter++;
          }
  }
  

/***  
  if(pid == e){
      felectroncounter++;
      //felectronmom += aTrack->GetDynamicParticle()->GetTotalMomentum();
  }
  
  if(pid == pos){
      fposcounter++;
      //felectronmom += aTrack->GetDynamicParticle()->GetTotalMomentum();
  }
  if(pid == y){
      fgammacounter++;
      //fgammamom += aTrack->GetDynamicParticle()->GetTotalMomentum();
  }
  if(pid == antinu_mu){
     fnumucounter++;
     //fnumumom += aTrack->GetDynamicParticle()->GetTotalMomentum();
  }
  if(pid == nu_e){
     fantinuecounter++;
     //fantinuemom += aTrack->GetDynamicParticle()->GetTotalMomentum();
  }
  //else{
  //	G4cout << "Other particle: " << pid->GetParticleName() << G4endl;
  //}

***/

  else{
//  if(pid != pos && pid != y && pid != antinu_mu && pid != nu_e  && pid != opticalphoton && pid != e){
  if(aTrack->GetTrackID() > 1) {
  if( (aTrack->GetGlobalTime() / CLHEP::ns) >= 3){ 
  G4cout << "Track ID: " << aTrack->GetTrackID() << G4endl;
  G4cout << "--------> PID: " << pid->GetParticleName() << G4endl;
  //G4cout << "--------> Process: " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
  //G4cout << "Momentum: " << aTrack->GetDynamicParticle()->GetTotalMomentum()  << "MeV" << G4endl;
  G4cout << "--------> Global Time: " << aTrack->GetGlobalTime() / CLHEP::ns << " ns" << G4endl;
  // Position
  G4ThreeVector pos = aTrack->GetPosition();
  G4double z = pos.getZ();

  // Origin
  G4ThreeVector vertex = aTrack->GetVertexPosition();
  G4double origin = vertex.getZ();

  G4double tracklength = z - origin;
  G4cout << "Position Z : " << tracklength/ CLHEP::mm << " mm" << G4endl;
  //}
  }
  }
  }
  return fUrgent;
}

//....

void StackingAction::NewStage()
{


  auto analysisManager = G4AnalysisManager::Instance();

  if(fCerenkovCounter>0) {
  analysisManager->FillH1(3, fCerenkovCounter);
  }

  //if(fScintillationCounter>0) {
  //analysisManager->FillH1(8, fScintillationCounter);
  //}

  G4cout << "Number of Scintillation photons produced in this event : "
         << fScintillationCounter << G4endl;
  G4cout << "Number of Cerenkov photons : "
	 << fCerenkovCounter << G4endl;

/***  
  G4cout << "------> Electrons <--------" << G4endl;
  G4cout << "Number: " << felectroncounter << G4endl;
  //G4cout << "Total momentum: " << felectronmom << G4endl;

  G4cout << "------> Muon Neutrinos <-------" << G4endl;
  G4cout << "Number: " << fnumucounter << G4endl;
  //G4cout << "Total momentum: " << fnumumom << G4endl;

  G4cout << "------> Anti-electron Neutrinos------" << G4endl;
  G4cout << "Number: " << fantinuecounter << G4endl;
  //G4cout << "Total momentum: " << fantinuemom << G4endl;

  G4cout << "------> Gammas <-------" << G4endl;
  G4cout << "Number: " << fgammacounter << G4endl;
  //G4cout << "Total momentum: " << fgammamom << G4endl;
***/

}

//....

void StackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  felectroncounter = 0;
  fposcounter = 0;
  fgammacounter = 0;
  fnumucounter = 0;
  fantinuecounter = 0;

}

//....

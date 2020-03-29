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

#include "TrackingAction.hh"

//#include "Run.hh"
#include "EventAction.hh"
//#include "TrackingMessenger.hh"
#include "Analysis.hh"
#include "RunAction.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//.....

TrackingAction::TrackingAction()
:G4UserTrackingAction()
 
{
  
  fTimeWindow1 = fTimeWindow2 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::SetTimeWindow(G4double t1, G4double dt)
{
 // fTimeWindow1 = t1;
  //fTimeWindow2 = fTimeWindow1 + dt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  //Run* run 
  // = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
         
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String name   = particle->GetParticleName();
    
  //count particles
  //run->ParticleCount(name, Ekin, meanLife);
  
  fTime_birth = track->GetGlobalTime();
}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  
 //Run* run 
 // = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
   
 // When Primary Particle Decays 
 G4int ID = track->GetTrackID();
 if (ID==1 && track->GetKineticEnergy() == 0) { 
     G4double time = track->GetGlobalTime() / ns; 
     G4cout << "TIME " << time << "ns" << G4endl;

     // Position
     G4ThreeVector pos = track->GetPosition();
     G4double z = pos.getZ();

     // Origin
     G4ThreeVector vertex = track->GetVertexPosition();
     G4double origin = vertex.getZ();

     G4double tracklength = z - origin;
 
     G4cout << "Delta z : " << tracklength/mm << " mm" << G4endl;

     //G4cout << "POSITION DECAY " << z/mm << G4endl;
     G4AnalysisManager* analysis = G4AnalysisManager::Instance();

     analysis->FillH1(7, time); 
     analysis->FillH1(8, tracklength/mm);

 }

 //  G4int ID = track->GetTrackID();
 // if (ID == 1) run->PrimaryTiming(time);        //time of life of primary ion
 // fTime_end = time;      

 
 //count activity in time window
 //run->SetTimeWindow(fTimeWindow1, fTimeWindow2);
  
 //G4bool life1(false), life2(false), decay(false);
 //if ((fTime_birth <= fTimeWindow1)&&(fTime_end > fTimeWindow1)) life1 = true;
 //if ((fTime_birth <= fTimeWindow2)&&(fTime_end > fTimeWindow2)) life2 = true;
 //if ((fTime_end   >  fTimeWindow1)&&(fTime_end < fTimeWindow2)) decay = true;
 //if (life1||life2||decay) run->CountInTimeWindow(name,life1,life2,decay);
}

//....

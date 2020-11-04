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
// nnbar calorimeter prototype simulation for GEANT4 tutorial 
//
// Author: Katherine Dunne
// email: katherine.dunne@fysik.su.se
//
///////////////////////////////////////////////////////////////////////


#include "EventAction.hh"
#include "Analysis.hh"
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "NNbarHit.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....

EventAction::EventAction(): 
    G4UserEventAction(),
    fAbsoEdepHCID(-1),
    fGapEdepHCID(-1),
    fAbsoTrackLengthHCID(-1),
    fCerenkovHCID(-1),
    fScintTrackLengthHCID(-1),
    scintHitsCollectionID(-1),
    absHitsCollectionID(-1)
{}

//....

EventAction::~EventAction()
{}

//....

G4THitsMap<G4double>* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<G4THitsMap<G4double>*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....

G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;  
}  

//.....

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(scintHitsCollectionID == -1) {
       scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
  }


}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  


   G4HCofThisEvent* HCE = event->GetHCofThisEvent();
 
    if(scintHitsCollectionID  < 0) {
        return;
    }

        int CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");
        int CHCID2 = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");

    NNbarHitsCollection* ScintHits = 0;
    NNbarHitsCollection* AbsHits   = 0;

    if (HCE) {
        G4AnalysisManager* analysis = G4AnalysisManager::Instance();
	G4int ltime     = 0.;
    	G4int parentID  = 0;
	G4String proc   = "";
	G4String name   = "";
        G4double time   = 0.;
        G4int trID      = 0;
	G4int i         = 0;
	G4double kinEn  = 0.;
	G4double eDep   = 0.;
        G4double trackl = 0.;	
        G4int hitCount  = 0;

        // Book vector to keep track of Edep in each Scintillator Sheet
	G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
        G4double totEdep   = 0.;     
        G4double eDepScint = 0.;
        G4double eDepAbs   = 0.;
        G4int cerenkovCounter = 0;

	ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID));
	hitCount = ScintHits->entries();

        if (hitCount > 0) {
	    for (G4int h=0; h<hitCount; h++) {
	        ltime    = ((*ScintHits)[h]) -> GetLocalTime();
	        parentID = ((*ScintHits)[h]) -> GetParentID();
    		proc     = ((*ScintHits)[h]) -> GetProcess();
       	        name     = ((*ScintHits)[h]) -> GetName();
       	        time     = ((*ScintHits)[h]) -> GetTime(); 
	        trID     = ((*ScintHits)[h]) -> GetTrackID();
		i        = ((*ScintHits)[h]) -> GetXID();
	        kinEn    = ((*ScintHits)[h]) -> GetKinEn();
	        eDep     = ((*ScintHits)[h]) -> GetEdep();
                trackl   = ((*ScintHits)[h]) -> GetPosZ();	

               if (proc != "primary") {
                   continue;
               }

               // Sum totEdep
               eDepScint += eDep;  
               totEdep += eDep;
                

	    }
	       
            // Fill total Edep in Scintillator
            if (eDepScint > 0) {
                analysis->FillH1(14, eDepScint/CLHEP::MeV);	
            }
      }

      AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID2));
      hitCount = AbsHits->entries();

      if(hitCount > 0) {
  
          for (G4int h=0; h<hitCount; h++) {
	      ltime     = ((*AbsHits)[h]) -> GetLocalTime();
	      parentID	= ((*AbsHits)[h]) -> GetParentID();
     	      proc      = ((*AbsHits)[h]) -> GetProcess();
	      name      = ((*AbsHits)[h]) -> GetName();
	      time      = ((*AbsHits)[h]) -> GetTime(); 
	      trID      = ((*AbsHits)[h]) -> GetTrackID();
	      i         = -99;
	      kinEn     = ((*AbsHits)[h]) -> GetKinEn();
	      eDep      = ((*AbsHits)[h]) -> GetEdep();
              trackl    = ((*AbsHits)[h]) -> GetPosZ();	
                
              if (proc != "primary") {
                    continue;
              }
               
              eDepAbs += eDep;
              totEdep += eDep;
                
              if (name == "opticalphoton" && proc == "Cerenkov" && ltime == 0){
	          analysis->FillH1(11,time);
	          cerenkovCounter++;
	      }
               
                              	    
	   }
           
           if (eDepAbs>0) {
               analysis->FillH1(15, eDepAbs/CLHEP::MeV);
           }
  
    }
       
 
    G4cout << "Energy Deposited in Scintillator: " << eDepScint/CLHEP::MeV << " MeV" << G4endl;
    G4cout << "Energy Deposited in Absorber: " << eDepAbs/CLHEP::MeV << " MeV" << G4endl; 
      
  } else {
        G4cout << "No HCE" << G4endl;
  }

  //print per event (modulo n)
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     
   // PrintEventStatistics(absoEdep, absoTrackLength, gapEdep, scintTrackLength);
  }
}  

//....

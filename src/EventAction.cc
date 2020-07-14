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


#include "EventAction.hh"
#include "Analysis.hh"
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "TubeSD.hh"
#include "NNbarHit.hh"
#include "TrackerSD.hh"

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
    absHitsCollectionID(-1),
    tubeHitsCollectionID(-1),
    trackerHitsCollectionID(-1)
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
       tubeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");  	
       trackerHitsCollectionID = pSDManager->GetCollectionID("TrackerHitsCollection"); 
 }


}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    if(scintHitsCollectionID  < 0) {
        return;
    }

    int CHCID = -1;
    if (CHCID<0) {
        CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");
    }
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();


    int CHCID2 = -1;
    if (CHCID2<0) {
        CHCID2 = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");
    }

    int CHCID3 = -1;
    if (CHCID3<0) {
        CHCID3 = G4SDManager::GetSDMpointer()->GetCollectionID("TubeHitCollection");
    }

    int CHCID4 = -1;
    if (CHCID4<0) {
	CHCID4 = G4SDManager::GetSDMpointer()->GetCollectionID("TrackerHitsCollection");
    }

    NNbarHitsCollection* ScintHits = 0;
    NNbarHitsCollection* AbsHits   = 0;
    NNbarHitsCollection* TubeHits  = 0;
    NNbarHitsCollection* TrackerHits = 0;

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

	ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID));

        // Book vector to keep track of Edep in each Scintillator Sheet
	G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
        G4double totEdep   = 0.;     
        G4double eDepScint = 0.;
        G4double eDepAbs   = 0.;
        G4double eDepTube  = 0.; 
	G4double eDepTracker = 0.;

        if (ScintHits) {
	    hitCount = ScintHits->entries();

	    for (G4int h=0; h<hitCount; h++) {
	        // In future can instead define aHit and assign like track.member[i]
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

               // Sum totEdep
               eDepScint += eDep;  
               totEdep += eDep;
               if (trID ==1) {
                  // Sum eDep for each scintillator sheet
                   EdepPerSheet[i] += eDep;
	
		   analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	           if (kinEn == 0) {
                       if(h==0){
                           G4cout << "Filling with pos " << trackl << G4endl;
                           analysis->FillH1(13, trackl/CLHEP::cm);
		           analysis->FillH1(12, time/CLHEP::ns);
                           continue;
                       }
		       
                       //G4cout << "hit number: " << h << G4endl;
                       G4double prevKin = ((*ScintHits)[h-1])->GetKinEn();
                       //G4cout << "Position: " << trackl << G4endl;
	               //G4cout << "Previous KinEn: " << prevKin << G4endl;
		       //G4cout << "Local Time: " << ltime << G4endl;
		       
                       // For some reason, kinEn can be 0 two hits in a row-> double counting one primary particle
                       if (prevKin == 0) {
                           G4cout << "continuing to next hit "<< G4endl;
                           continue;
                       }
 
			else {
                           //G4cout << "Filling with pos " << trackl << G4endl;
                           analysis->FillH1(13, trackl/CLHEP::cm);
		           analysis->FillH1(12, time/CLHEP::ns);
                           analysis->FillH2(1, trackl/CLHEP::cm, eDepScint/CLHEP::MeV);
                       }
		   }
	       } 

	    }

	    // Fill scint bins with Energy Dep
            for (G4int i=0; i<10; i++) {
                if (EdepPerSheet[i]) {
		    analysis->FillH1(i, EdepPerSheet[i]/CLHEP::MeV);
                }
	    }
   
            // Fill total Edep in Scintillator
            analysis->FillH1(14, eDepScint/CLHEP::MeV);	
            G4cout << "Total Edep in scint: " << eDepScint/CLHEP::MeV << G4endl;         
        }
 
	AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID2));
       

	if(AbsHits) {
	
	    hitCount = AbsHits->entries();
            G4int cerenkovCounter = 0;
	    for (G4int h=0; h<hitCount; h++) {
		ltime           = ((*AbsHits)[h]) -> GetLocalTime();
		parentID	= ((*AbsHits)[h]) -> GetParentID();
     	        proc            = ((*AbsHits)[h]) -> GetProcess();
	        G4String name   = ((*AbsHits)[h]) -> GetName();
	    	G4double time   = ((*AbsHits)[h]) -> GetTime(); 
		G4int trID      = ((*AbsHits)[h]) -> GetTrackID();
		G4int i         = -99;
	        G4double kinEn  = ((*AbsHits)[h]) -> GetKinEn();
	        G4double eDep   = ((*AbsHits)[h]) -> GetEdep();
                G4double trackl = ((*AbsHits)[h]) -> GetPosZ();	
                
 
                //if(trID>1 && eDep>0){
                    //G4cout << "Particle: " << name << "   Process: "<< proc << "   Deposit: " << eDep/CLHEP::MeV << G4endl;
                //}
                if (proc == "Decay") {
                    continue;
                }
                	   
                eDepAbs += eDep;
                totEdep += eDep;
                if (trID == 1){
                    analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	            if (kinEn == 0) {
                        if(h==0){
                            analysis->FillH1(13, trackl/CLHEP::cm);
		            analysis->FillH1(12, time/CLHEP::ns);
                            continue;
                        }
                        G4double prevKin = ((*AbsHits)[h-1])->GetKinEn();
                        if (prevKin==0){
                            //G4cout << "continuing to next hit" << G4endl;
                            continue;
                        } else {
                            analysis->FillH1(13, trackl/CLHEP::cm);
		            analysis->FillH1(12, time/CLHEP::ns);
                            analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);
                       }
		    }
	        }
                // ltime == 0 ?
                if (name == "opticalphoton" && proc == "Cerenkov"){
                    //G4cout << "photon local time: " << ltime << G4endl;
		    analysis->FillH1(11,time);
		    cerenkovCounter++;
		}
	    
	    }
	    
            if (eDepAbs>0){
                analysis->FillH1(15, eDepAbs/CLHEP::MeV);
                //analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);  
                  
                if(cerenkovCounter>0){
                    G4cout << "Filling ceren histos with " << cerenkovCounter << " photons and " << eDepAbs/CLHEP::MeV << " MeV deposited." << G4endl;
                    analysis->FillH2(2, cerenkovCounter, eDepAbs/CLHEP::MeV);       
                }
            }

            G4cout << "Total Edep in lead-glass: " << eDepAbs/CLHEP::MeV << G4endl;
	    G4cout << "Cerenkov count: " << cerenkovCounter << G4endl;
            if (cerenkovCounter>0){
	        analysis->FillH1(10,cerenkovCounter);
	    }
	}  
       
        TrackerHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID4));
 
        if (TrackerHits) {
	    hitCount = TrackerHits->entries();
            G4cout << " in trackerhits loops " << G4endl;
	    for (G4int h=0; h<hitCount; h++) {
//	        // In future can instead define aHit and assign like track.member[i]
	        ltime    = ((*TrackerHits)[h]) -> GetLocalTime();
//	        parentID = ((*TrackerHits)[h]) -> GetParentID();
//	        proc     = ((*TrackerHits)[h]) -> GetProcess();
//              name     = ((*TrackerHits)[h]) -> GetName();
//              time     = ((*TrackerHits)[h]) -> GetTime(); 
//	        trID     = ((*TrackerHits)[h]) -> GetTrackID();
//		i        = ((*TrackerHits)[h]) -> GetXID();
//	        kinEn    = ((*TrackerHits)[h]) -> GetKinEn();
//	        eDep     = ((*TrackerHits)[h]) -> GetEdep();
//              trackl   = ((*TrackerHits)[h]) -> GetPosZ();	
 
               // Sum totEdep
//               eDepTracker += eDep;  
//               totEdep += eDep;
//               if (trID ==1) {
//		   analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
//	           if (kinEn == 0) {
//                       if(h==0){
//                           G4cout << "Filling with pos " << trackl << G4endl;
//                           analysis->FillH1(13, trackl/CLHEP::cm);
//		           analysis->FillH1(12, time/CLHEP::ns);
//                           continue;
//                       }
		       
//                       //G4cout << "hit number: " << h << G4endl;
//                       G4double prevKin = ((*TrackerHits)[h-1])->GetKinEn();
//                       //G4cout << "Position: " << trackl << G4endl;
//	               //G4cout << "Previous KinEn: " << prevKin << G4endl;
//		       //G4cout << "Local Time: " << ltime << G4endl;
		       
//                       // For some reason, kinEn can be 0 two hits in a row-> double counting one primary particle
//                       if (prevKin == 0) {
//                           G4cout << "continuing to next hit "<< G4endl;
//                           continue;
//                       } else {
//                           //G4cout << "Filling with pos " << trackl << G4endl;
//                           analysis->FillH1(13, trackl/CLHEP::cm);
//		           analysis->FillH1(12, time/CLHEP::ns);
//                           //analysis->FillH2(1, trackl/CLHEP::cm, eDepScint/CLHEP::MeV);
//                       }
//		   }
//	       } 
//
	    }
   
//            // Fill total Edep in Trackere
//            analysis->FillH1(16, eDepTracker/CLHEP::MeV);	
//            G4cout << "Total Edep in tracker: " << eDepTracker/CLHEP::MeV << G4endl;         
        }

        TubeHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID3));
 
        if (TubeHits) {
	    hitCount = TubeHits->entries();
            G4cout << " in tubehits loops " << G4endl;
	    for (G4int h=0; h<hitCount; h++) {
	        ltime    = ((*TubeHits)[h]) -> GetLocalTime();
	        parentID = ((*TubeHits)[h]) -> GetParentID();
    		proc     = ((*TubeHits)[h]) -> GetProcess();
       	        name     = ((*TubeHits)[h]) -> GetName();
       	        time     = ((*TubeHits)[h]) -> GetTime(); 
	        trID     = ((*TubeHits)[h]) -> GetTrackID();
		i        = ((*TubeHits)[h]) -> GetXID();
	        kinEn    = ((*TubeHits)[h]) -> GetKinEn();
	        eDep     = ((*TubeHits)[h]) -> GetEdep();
                trackl   = ((*TubeHits)[h]) -> GetPosZ();

	        eDepTube += eDep;  
               totEdep += eDep;
               if (trID ==1) {
		   analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	           if (kinEn == 0) {
                       if(h==0){
                           G4cout << "Filling with pos " << trackl << G4endl;
                           analysis->FillH1(13, trackl/CLHEP::cm);
		           analysis->FillH1(12, time/CLHEP::ns);
                           continue;
                       }
                       G4double prevKin = ((*TubeHits)[h-1])->GetKinEn();

                       if (prevKin == 0) {
                           G4cout << "continuing to next hit "<< G4endl;
                           continue;
                       } else {

                           analysis->FillH1(13, trackl/CLHEP::cm);
		           analysis->FillH1(12, time/CLHEP::ns);
                       }
		   }
	       } 

	    }
   
            analysis->FillH1(16, eDepTube/CLHEP::MeV);	
            G4cout << "Total Edep in tube: " << eDepTube/CLHEP::MeV << G4endl;         
        }




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

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
#include "NNbarHit.hh"
//#include "AbsorberHit.hh"
//#include "ScintillatorHit.hh"


#include "G4VHitsCollection.hh"
//#include "G4VEventManager.hh"
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

/***void EventAction::PrintEventStatistics(
                            G4double absoEdep, G4double absoTrackLength,
                            G4double gapEdep, G4double scintTrackLength) const //, G4double numCerenkov) const
{
  // Print event statistics
  //
  //
  
  G4cout
     << "----->Absorber<-----" << G4endl  << "total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       Total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl;
  
  G4cout 
     << "------>Scint<-------" << G4endl << "total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       Total track length: " 
     << std::setw(7) << G4BestUnit(scintTrackLength, "Length")
     << G4endl;i

    
}
***/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(scintHitsCollectionID == -1) {
       scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  
   // Get hist collections IDs
  //if ( fAbsoEdepHCID == -1 ) {
    //fAbsoEdepHCID 
    //  = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber/Edep");
    //fGapEdepHCID 
    //  = G4SDManager::GetSDMpointer()->GetCollectionID("Scint/Edep");
   // fAbsoTrackLengthHCID 
    //  = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber/TrackLength");
    //fCerenkovHCID
    //  = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber/Population");
    //fScintTrackLengthHCID
    //  = G4SDManager::GetSDMpointer()->GetCollectionID("Scint/TrackLength");


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

    NNbarHitsCollection* ScintHits = 0;
    NNbarHitsCollection* AbsHits = 0;

    if (HCE) {
       //G4cout << "In HCE loop: " << G4endl;

	ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID));
        //ScintillatorHitsCollection* sCHC = (ScintillatorHitsCollection*)(HCE->GetHC(scintHitsCollectionID));
        //AbsorberHitsCollection* aCHC =(AbsorberHitsCollection*)(HCE->GetHC(absHitsCollectionID));

	G4double totEdep[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};

        if (ScintHits) {
	//if (sCHC) {
	    G4int hitCount = ScintHits->entries();
            G4AnalysisManager* analysis = G4AnalysisManager::Instance();
	    

            G4int trID   = 0;
	    G4int i         = 0;
	    G4double kinEn  = 0.;
	    G4double eDep   = 0.;
            G4double trackl = 0.;	
            

	    for (G4int h=0; h<hitCount; h++) {
	        // In future can instead define aHit and assign like track.member[i]
	        G4int trID   = ((*ScintHits)[h]) -> GetTrackID();
		G4int i         = ((*ScintHits)[h]) -> GetXID();
	        G4double kinEn  = ((*ScintHits)[h]) -> GetKinEn();
	        G4double eDep   = ((*ScintHits)[h]) -> GetEdep();
                G4double trackl = ((*ScintHits)[h]) -> GetPosZ();	
                
	        totEdep[i] += eDep;	
                   G4cout << "Scint Replica Number: " << i << G4endl;
	           G4cout << "Scint Kinetic Energy: " << kinEn << G4endl;
	           G4cout << "Scint Energy Deposited: " << eDep << G4endl;
	           G4cout << "Scint Position: " << trackl/CLHEP::cm << " cm" << G4endl;

                
		if (trID == 1) {
		    analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
                    if (kinEn == 0) {
		        analysis->FillH1(13, trackl/CLHEP::cm);
	            }
		}
	    }


	    // Fill scint bins with Energy Dep

            for (G4int i=0; i<10; i++) {
                if (totEdep[i]) {
		    G4cout << "totEdep in Scint " << i <<  " : " << totEdep[i]/CLHEP::MeV << G4endl;  
		    analysis->FillH1(i, totEdep[i]/CLHEP::keV);
                }
	    }
          }
	

	AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID2));

	if(AbsHits) {
	
	    G4int hitCount = AbsHits->entries();

            G4AnalysisManager* analysis = G4AnalysisManager::Instance();
	    
	    for (G4int h=0; h<hitCount; h++) {
	        G4int trID      = ((*AbsHits)[h]) -> GetTrackID();
		G4int i         = -99;
	        G4double kinEn  = ((*AbsHits)[h]) -> GetKinEn();
	        G4double eDep   = ((*AbsHits)[h]) -> GetEdep();
                G4double trackl = ((*AbsHits)[h]) -> GetPosZ();	

	        G4cout << "Absorber Kinetic Energy: " << kinEn << G4endl;
	        G4cout << "Absorber Energy Deposited: " << eDep << G4endl;
	        G4cout << "Absorber Position: " << trackl/CLHEP::cm << " cm" << G4endl;
                

	        if (trID ==1){
                    analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	            if (kinEn == 0) {
		        analysis->FillH1(13, trackl/CLHEP::cm);
		    }
	        }
	    
	    }
	}  
        else {
	    G4cout << "No Hits" << G4endl;
	}
    } else {
        G4cout << "No HCE" << G4endl;
    }


     
  //}
  
  // Get sum values from hits collections
/***
auto absoEdep = GetSum(GetHitsCollection(fAbsoEdepHCID, event));
  auto gapEdep = GetSum(GetHitsCollection(fGapEdepHCID, event));
  auto absoTrackLength 
    = GetSum(GetHitsCollection(fAbsoTrackLengthHCID, event));
  auto numcerenkov = GetSum(GetHitsCollection(fCerenkovHCID, event));
  auto scintTrackLength = GetSum(GetHitsCollection(fScintTrackLengthHCID, event));
  auto totTrackLength = scintTrackLength + absoTrackLength;
***/

  // get analysis manager
  //auto analysisManager = G4AnalysisManager::Instance();

 
  // fill histograms
  



/***

  if (absoEdep >0 ) {
  analysisManager->FillH1(0, absoEdep);
  }

  if (gapEdep >0) {
  analysisManager->FillH1(1, gapEdep);
  }
  if (absoTrackLength >0) {
  analysisManager->FillH1(2, absoTrackLength);
  }
  if(numcerenkov > 0) {
      analysisManager->FillH1(3, numcerenkov);
  }
  ***/

  //analysisManager->FillH1(4, scintTrackLength);
  //analysisManager->FillH1(5, totTrackLength);
 


  //print per event (modulo n)
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     
   // PrintEventStatistics(absoEdep, absoTrackLength, gapEdep, scintTrackLength);
  }
}  

//....

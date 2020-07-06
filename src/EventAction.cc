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
    tubeHitsCollectionID(-1)
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
         tubeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");  
         scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
         absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");

 }


}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    if(scintHitsCollectionID  < 0) {
        return;
    }

    G4HCofThisEvent* HCE = event->GetHCofThisEvent();


    int CHCID = -1;
    if (CHCID<0) {
        CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");
    }


    int CHCID2 = -1;
    if (CHCID2<0) {
        CHCID2 = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");
    }

    int CHCID3 = -1;
    if (CHCID3<0) {
        CHCID3 = G4SDManager::GetSDMpointer()->GetCollectionID("TubeHitCollection");
    }


    //G4cout << "ScintHCID" << CHCID << G4endl;
    //G4cout << "TubeHCID" << CHCID3 << G4endl;
    //G4cout << "AbsorberHCID" << CHCID2 << G4endl;

    NNbarHitsCollection* ScintHits = 0;
    NNbarHitsCollection* AbsHits   = 0;
    NNbarHitsCollection* TubeHits  = 0;
    NNbarHitsCollection* Hits = 0;


    // Initialize to -99 so can see if real 0 or not
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    G4double ltime     = -99.;
    G4int parentID  = -99;
    G4String proc   = "";
    G4String name   = "";
    G4double time   = -99.;
    G4int trID      = -99;
    G4int vox       = -99;
    G4double kinEn  = -99.;
    G4double eDep   = -99.;
    G4double trackl = -99.;	
    G4int hitCount  = -99;

    // Book vector to keep track of Edep in each Scintillator Sheet
    G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
    G4double totEdep          = -99.;     
    G4double eDepScint        = -99.;
    G4double eDepAbs          = -99.;
    G4double eDepTube         = -99.; 
    G4double eDepElectron     = -99.;
    G4double primaryStop      = -99.;
    G4double eDepInel         = -99.;
    G4double eDepIoni         = -99.;
    G4double eDepComp         = -99.;
    G4double eDepPrimary      = -99.;
    G4int cerenkovCounter     = -99;


    if (HCE) {

        G4String HCName[3] = {"Tube", "Scint", "Absorber"};

        for (int i=0; i<3; i++) {
            // Loop through each hit collection
            Hits = (NNbarHitsCollection*)(HCE->GetHC(i));
            hitCount = Hits->entries();

            G4cout << "Looping through " << hitCount << " " << HCName[i] << " Hits." << G4endl;

	    if (i == 2 && hitCount > 0) {
                cerenkovCounter = 0;
            }
            G4cout << "numCerenkov " << cerenkovCounter << G4endl;

            for (G4int h=0; h<hitCount; h++) {
	        // In future can instead define aHit and assign like track.member[i]
	        ltime    = ((*Hits)[h]) -> GetLocalTime();
	        parentID = ((*Hits)[h]) -> GetParentID();
    		proc     = ((*Hits)[h]) -> GetProcess();
       	        name     = ((*Hits)[h]) -> GetName();
       	        time     = ((*Hits)[h]) -> GetTime(); 
	        trID     = ((*Hits)[h]) -> GetTrackID();
		vox        = ((*Hits)[h]) -> GetXID();
	        kinEn    = ((*Hits)[h]) -> GetKinEn();
	        eDep     = ((*Hits)[h]) -> GetEdep();
                trackl   = ((*Hits)[h]) -> GetPosZ();	
                //G4cout << "Name : " << name << G4endl;

                if (proc != "primary" && proc !="Cerenkov" || time > 10*CLHEP::ns){
                    continue;
                }

                //if (proc != "primary") continue;


                //if (proc == "Decay"){
                //    continue;
                //}

                // Sum totEdep
                totEdep += eDep;

                if(proc == "protonInelastic") {
                  eDepInel += eDep;
                } else if (proc == "hIoni"){
                  eDepIoni += eDep;
                } else if (proc == "compt") {
                  eDepComp += eDep;
                } else {
                  eDepPrimary += eDep;
                }

 
                // Record tracklength when particle and stop tracking when stops
                if (trID ==1) {
                   analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);


	           if (kinEn == 0) {
                       //if(h==0){
                           primaryStop = trackl;
                           
                       //    continue;
                       //}
		       
                       //G4cout << "hit number: " << h << G4endl;
                       //G4double prevKin = ((*Hits)[h-1])->GetKinEn();
		       
                       // For some reason, kinEn can be 0 two hits in a row-> double counting one primary particle
                       //if (prevKin == 0) {
                           //G4cout << "continuing to next hit "<< G4endl;
                       //    continue;
                       //} else {
                       //    primaryStop = trackl;
                       //}
		   }

	        } // End primary stop loop

                
                if (i==0) { // Vacuum Tube Hits
                    eDepTube += eDep;  

                } else if (i==1) {  // Scintillator Hits
                    EdepPerSheet[vox] += eDep;
                    eDepScint += eDep;
                } else if (i==2) { // Absorber Hits
                    eDepAbs += eDep;
                    
                    // ltime==0 don't double count same photon
                    if (ltime == 0 && name == "opticalphoton" && proc == "Cerenkov"){
                        if (primaryStop/CLHEP::cm < 32 && primaryStop/CLHEP::cm > -9.9) {
                            //G4cout << "time: " << time/CLHEP::ns << G4endl;
                           // if (time/CLHEP::ns > 100) {
                            //    G4cout << "name: " << name << G4endl;
                           //     G4cout << "proc: " << proc << G4endl;
                           //     G4cout << "trackl: " << trackl/CLHEP::cm << G4endl;
                           // }
                            //G4cout << "trackl: " << trackl/CLHEP::cm << G4endl;
                            //G4cout << "primaryStop: " << primaryStop/CLHEP::cm << G4endl;
                        }
		        analysis->FillH1(11,time);
		        cerenkovCounter++;
                    }
	
                } else {
                    G4cout << "incorrect HCID" << G4endl;
                }


            } // End hitCount loop

            //G4cout << " end of loop" << G4endl;

        } // End HCName loop
   


        /////////////////////////
        //
        // Fill Histograms
        //
        /////////////////////////

 
        // Fill total Edep in Vacuum Tube
        if (eDepTube > 0) {
            analysis->FillH1(16, eDepTube/CLHEP::MeV);	
        } 

        // Fill scint bins with Energy Dep
        for (G4int j=0; j<10; j++) {
            if (EdepPerSheet[j]) {
	        analysis->FillH1(j, EdepPerSheet[j]/CLHEP::MeV);
            }
        }
   
        // Fill total Edep in Scintillator
        if (eDepScint>0) {
            analysis->FillH1(14, eDepScint/CLHEP::MeV);	
        }            

        //if (eDepAbs>0) {
            analysis->FillH1(15, eDepAbs/CLHEP::MeV);
            //analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);  
                  
            if(eDepAbs>-99 && cerenkovCounter>-99){
                             // G4cout << "Filling ceren histos with " << cerenkovCounter << " photons and " << eDepAbs/CLHEP::MeV << " MeV deposited." << G4endl;
                
                analysis->FillH2(5, cerenkovCounter, primaryStop/CLHEP::cm);
                analysis->FillH1(10, cerenkovCounter);                
            }
        //}

        analysis->FillH1(13, primaryStop/CLHEP::cm); 
        analysis->FillH2(1, primaryStop/CLHEP::cm, totEdep/CLHEP::MeV);

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

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
       scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
       tubeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");  
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

    NNbarHitsCollection* ScintHits = 0;
    NNbarHitsCollection* AbsHits   = 0;
    NNbarHitsCollection* TubeHits  = 0;

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
        G4double posZ   = 0.;
        G4double posX   = 0.;
        G4double posY   = 0.;
        G4double vertZ   = 0.;
        G4double vertX   = 0.;
        G4double vertY   = 0.;


        G4double trackl = 0.;	
        G4int hitCount  = 0;

	ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID));

        // Book vector to keep track of Edep in each Scintillator Sheet
	G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
        G4double totEdep   = 0.;     
        G4double eDepScint = 0.;
        G4double eDepAbs   = 0.;
        G4double eDepTube  = 0.; 
        G4double extraEdep = 0.;
        G4double eDepCompt = 0.;
        G4double eDepInelastic= 0.;
        G4double eDephIoni = 0.;
        G4double eDepHadElas = 0.;
        G4double eDepPrimary = 0.;
        G4double eDepOther = 0.;
        G4bool stopped = 0;
        G4bool de_stopped = 0;
        G4int trackID = 0;

        if (ScintHits) {
	    hitCount = ScintHits->entries();
            std::cout << "Processing: " << hitCount << " Scintillator Hits" << std::endl;
	    for (G4int h=0; h<hitCount; h++) {

                //if (h==0) G4cout << "Hit eDep: " << eDep << G4endl;
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
                posX     = ((*ScintHits)[h]) -> GetPosX();	
                posY     = ((*ScintHits)[h]) -> GetPosY();	
                posZ     = ((*ScintHits)[h]) -> GetPosZ();	
                vertX     = ((*ScintHits)[h]) -> GetVertX();	
                vertY     = ((*ScintHits)[h]) -> GetVertY();	
                vertZ     = ((*ScintHits)[h]) -> GetVertZ();	


	      if (name=="deuteron") {
                   if (kinEn == 0 && de_stopped!=1){
                       std::cout << " Stopping position: " << posZ/CLHEP::cm << std::endl;
                       de_stopped=1;
                   } 
                   //std::cout << "deuteron KE: " << kinEn/CLHEP::MeV << std::endl;
                   //std::cout << "posZ: " << posZ << std::endl;
		   //std::cout << "proc: " << proc << std::endl;
                   //std::cout << "trackID: " << trID << std::endl;
                   if (trackID!=trID) {
                   analysis->FillH1(17, kinEn/CLHEP::MeV);
                   analysis->FillH1(18, vertZ/CLHEP::cm);
                   trackID = trID;
                   }
              }


                //if(trID>1 && eDep>0){
                  //G4cout << "Particle: " << name << "   Process: "<< proc << "   Deposit: " << eDep/CLHEP::MeV << G4endl;
               //}
               if (proc == "Decay") {
                   continue;
               }


               // Sum totEdep
               eDepScint += eDep;  
               totEdep += eDep;
/***                

               if (proc != "primary" & eDep > 0) {

                   G4cout << " Scintillator Edep by non primary particle" << G4endl;
                   G4cout << "----------------------------" << G4endl;
                   G4cout << "Particle: " << name << G4endl;
                   G4cout << "Edep: " << eDep/CLHEP::MeV << " MeV" << G4endl;
                   G4cout << "Process: " << proc << G4endl;
                   G4cout << "time: " << time/CLHEP::ns << " ns" << G4endl << G4endl;
		   extraEdep += eDep;
                   if (proc == "compt") eDepCompt += eDep;
		   else if (proc == "pi+Inelastic") eDepInelastic += eDep;
                   else if (proc == "hIoni") eDephIoni += eDep;
                   else if (proc == "hadElastic") eDepHadElas += eDep;
                   else eDepOther += eDep;
                   continue;
               }
***/

                  
               //G4cout << "eDepScint: " << eDepScint << G4endl;
               //G4cout << "totEdep: " << totEdep << G4endl; 

               // Sum eDep for each scintillator sheet
               EdepPerSheet[i] += eDep;

               if (trID ==1) {
                   // Sum eDep for each scintillator sheet
                   EdepPerSheet[i] += eDep;
	           eDepPrimary += eDep;

                   G4cout << "Kinetic Energy: " << kinEn << " " << G4endl;
                   G4cout << "Position Z: " << posZ/CLHEP::cm << std::endl;

		   analysis->FillH2(0, posZ/CLHEP::cm - 166., kinEn/CLHEP::MeV);
	           
                   // When primary particle stops
                   if (kinEn == 0 && stopped!=1) {
                       stopped = 1;
                       G4cout << "Filling with position: " << posZ/CLHEP::cm - 166. << " cm" << G4endl;
                       //G4cout << "Filling with time: " << time/CLHEP::ns << " ns" <<  G4endl;
                       //G4cout << "Filling pos vs eDepScint with energy: " << eDepScint/CLHEP::MeV << " MeV" << G4endl;
                       analysis->FillH1(13, posZ/CLHEP::cm - 166.);
		       analysis->FillH1(12, time/CLHEP::ns);
                       // Filling only when trackID==1
                       analysis->FillH2(1, posZ/CLHEP::cm, totEdep/CLHEP::MeV);
		       
                       //G4cout << "hit number: " << h << G4endl;
                       //G4double prevKin = ((*ScintHits)[h-1])->GetKinEn();
                       //G4cout << "Position: " << trackl << G4endl;
	               //G4cout << "Previous KinEn: " << prevKin << G4endl;
		       //G4cout << "Local Time: " << ltime << G4endl;
		       
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

                        
//            for (G4int j=0; j<10; j++) {
//                eDepPrimary += EdepPerSheet[j];
//            }


/***
           // if (eDepScint > 240){
                G4cout << "Total Edep in scint: " << eDepScint/CLHEP::MeV << G4endl;         
                G4cout << "Sheet 1: " << EdepPerSheet[0]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 2: " << EdepPerSheet[1]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 3: " << EdepPerSheet[2]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 4: " << EdepPerSheet[3]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 5: " << EdepPerSheet[4]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 6: " << EdepPerSheet[5]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 7: " << EdepPerSheet[6]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 8: " << EdepPerSheet[7]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 9: " << EdepPerSheet[8]/CLHEP::MeV << G4endl;
                G4cout << "Sheet 10: " << EdepPerSheet[9]/CLHEP::MeV << G4endl;
***/
            //}

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
           /***
                if (proc != "primary" & eDep > 0) {
                   G4cout << "Absorber Edep by non primary particle" << G4endl;
                   G4cout << "----------------------------" << G4endl;
                   G4cout << "Particle: " << name << G4endl;
                   G4cout << "Edep: " << eDep/CLHEP::MeV << " MeV" << G4endl;
                   G4cout << "Process: " << proc << G4endl;
                   G4cout << "time: " << time/CLHEP::ns << " ns" << G4endl << G4endl;
		   extraEdep += eDep;
                   if (proc == "compt") eDepCompt += eDep;
		   else if (proc == "pi+Inelastic") eDepInelastic += eDep;
                   else if (proc == "hIoni") eDephIoni += eDep;
                   else if (proc == "hadElastic") eDepHadElas += eDep;
                   else eDepOther += eDep;  
                   continue;
               }
***/
 	   
               
                if (trID == 1){
                    eDepPrimary += eDep;
                    analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	            if (kinEn == 0) {
                        //analysis->FillH1(13, trackl/CLHEP::cm);
		        analysis->FillH1(12, time/CLHEP::ns);
                        analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);
                    }
	        }
                
                // ltime == 0 ?
                if (name == "opticalphoton" && proc == "Cerenkov"){
                    //G4cout << "photon local time: " << ltime << G4endl;
		    analysis->FillH1(11,time);
		    cerenkovCounter++;
                    //G4cout << "Cerenkov counter: " << cerenkovCounter << G4endl;
		}
	    
	    }
	    
            if (eDepAbs>0){
                analysis->FillH1(15, eDepAbs/CLHEP::MeV);
                //analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);  
                  
                if(cerenkovCounter>0){
                    //G4cout << "Filling ceren histos with " << cerenkovCounter << " photons and " << eDepAbs/CLHEP::MeV << " MeV deposited." << G4endl;
                    analysis->FillH2(2, cerenkovCounter, eDepAbs/CLHEP::MeV);       
                }
            }

            //G4cout << "Total Edep in lead-glass: " << eDepAbs/CLHEP::MeV << G4endl;
	    //G4cout << "Cerenkov count: " << cerenkovCounter << G4endl;
            if (cerenkovCounter>0){
	        analysis->FillH1(10,cerenkovCounter);
	    }
	}  
       
        TubeHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID3));
 
        if (TubeHits) {
	    hitCount = TubeHits->entries();
            // G4cout << " in tubehits loops " << G4endl;
	    for (G4int h=0; h<hitCount; h++) {
	        // In future can instead define aHit and assign like track.member[i]
	        ltime    = ((*TubeHits)[h]) -> GetLocalTime();
	        parentID = ((*TubeHits)[h]) -> GetParentID();
    		proc     = ((*TubeHits)[h]) -> GetProcess();
       	        name     = ((*TubeHits)[h]) -> GetName();
       	        time     = ((*TubeHits)[h]) -> GetTime(); 
	        trID     = ((*TubeHits)[h]) -> GetTrackID();
		i        = ((*TubeHits)[h]) -> GetXID();
	        kinEn    = ((*TubeHits)[h]) -> GetKinEn();
	        eDep     = ((*TubeHits)[h]) -> GetEdep();
                posX     = ((*TubeHits)[h]) -> GetPosX();	
                posY     = ((*TubeHits)[h]) -> GetPosY();	
                posZ     = ((*TubeHits)[h]) -> GetPosZ();	
                vertX     = ((*TubeHits)[h]) -> GetVertX();	
                vertY     = ((*TubeHits)[h]) -> GetVertY();	
                vertZ     = ((*TubeHits)[h]) -> GetVertZ();	



  	      if (name=="deuteron") {
                   if (kinEn == 0 && de_stopped!=1){
                       std::cout << " Stopping position: " << posZ/CLHEP::cm << std::endl;
                       de_stopped=1;
                   } 
                   //std::cout << "deuteron KE: " << kinEn/CLHEP::MeV << std::endl;
                   //std::cout << "posZ: " << posZ << std::endl;
		   //std::cout << "proc: " << proc << std::endl;
                   //std::cout << "trackID: " << trID << std::endl;
                   if (trackID!=trID) {
                   analysis->FillH1(17, kinEn/CLHEP::MeV);
                   analysis->FillH1(18, vertZ/CLHEP::cm);
                   trackID = trID;
                   }
              }

                //if(trID>1 && eDep>0){
                    //G4cout << "Particle: " << name << "   Process: "<< proc << "   Deposit: " << eDep/CLHEP::MeV << G4endl;
                //}
                if (proc == "Decay") {
                    continue;
                }

                // Sum totEdep
                eDepTube += eDep;  
                totEdep += eDep;

/***

                if (proc != "primary" & eDep > 0) {
                   G4cout << "Tube Edep by non primary particle" << G4endl;
                   G4cout << "----------------------------" << G4endl;
                   G4cout << "Particle: " << name << G4endl;
                   G4cout << "Edep: " << eDep/CLHEP::MeV << " MeV" << G4endl;
                   G4cout << "Process: " << proc << G4endl;
                   G4cout << "time: " << time/CLHEP::ns << " ns" << G4endl << G4endl;
		   extraEdep += eDep;
                   if (proc == "compt") eDepCompt += eDep;
		   else if (proc == "pi+Inelastic") eDepInelastic += eDep;
                   else if (proc == "hIoni") eDephIoni += eDep;
                   else if (proc == "hadElastic") eDepHadElas += eDep;
                   else eDepOther += eDep;
	           continue;
               }

***/
                         
                if (trID ==1) {
	            eDepPrimary += eDep;
     	            analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	            if (kinEn == 0) {
             //           G4cout << "Filling with pos " << trackl << G4endl;
                        //analysis->FillH1(13, trackl/CLHEP::cm);
		        analysis->FillH1(12, time/CLHEP::ns);
		        analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);
                    }
	        }
	    } 
   
            // Fill total Edep in Vacuum Tube
            analysis->FillH1(16, eDepTube/CLHEP::MeV);	
            //G4cout << "Total Edep in tube: " << eDepTube/CLHEP::MeV << G4endl;         
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

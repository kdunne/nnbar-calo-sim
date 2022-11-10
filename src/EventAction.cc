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
#include <math.h>

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
/***
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(scintHitsCollectionID == -1) {
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
  }
***/

}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();

    int CHCID = -1;

    //G4String name[] = {"A"};
    G4String name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

    G4String scorer = "/Pop";


//  Scintillator Bar Hits
    for(int i=0; i<10; i++) {

      
      scorer = "/Pop";
      std::string scorerName = "scint_" + name[i] + scorer;

      CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
      //std::cout << "Got collection ID " << std::endl;
      //std::cout << "Name: " << scorerName << ": CHCID = " << CHCID << std::endl;
      auto pop = GetSum(GetHitsCollection(CHCID, event)); 

      //scorer = "/eDep";
      //scorerName = "scint_" + name[i] + scorer;
      //CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
      //std::cout << "Name: " << scorerName << " CHCID: " << CHCID << std::endl;
      //auto eDep = GetSum(GetHitsCollection(CHCID,event));


      //scorer = "/Pop";
      //scorerName = "fiber_" + name[i] + scorer;
      //CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
      //std::cout << "Name: " << scorerName << " CHCID: " << CHCID << std::endl;
      //auto popFiber = GetSum(GetHitsCollection(CHCID, event)); 

      //scorerName = "sipm_" + name[i] + scorer;
      //CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
      //std::cout << "Name: " << scorerName << "CHCID: " << CHCID << std::endl;
      //auto popSipm = GetSum(GetHitsCollection(CHCID, event)); 
        

      G4AnalysisManager* analysis = G4AnalysisManager::Instance();
      analysis->FillH1(i, pop);

      //if (eDep >0) { 
        //analysis->FillH1(i+1, eDep/CLHEP::MeV);
        //analysis->FillH2(i, eDep/CLHEP::MeV, pop);
        //analysis->FillH2(i+10, pop, popFiber);
       //}

      //std::cout << "Num SiPM Photons:  " << popSipm << std::endl;
      //std::cout << "Num Fiber Photons: " << popFiber << std::endl;
      std::cout << "Num Scint Photons in " << scorerName << " : " << pop << std::endl;
      //std::cout << "eDep: " << eDep/CLHEP::MeV << std::endl;

    }


/***
//  Fiber Hits
    for(int i=0; i<10; i++) {
       
        //std::cout << "name: " << name[i] << j << scorer << std::endl;
        scorer = "/Pop";
        std::string scorerName = "fiber_" + name[i] + scorer;

        int CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
        //std::cout << "Name: " << scorerName << "CHCID: " << CHCID << std::endl;
        auto pop = GetSum(GetHitsCollection(CHCID, event)); 

        scorer = "/eDep";
        scorerName = "fiber_" + name[i] + scorer;

        CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
        //std::cout << "Name: " << scorerName << "CHCID: " << CHCID << std::endl;
        auto eDep = GetSum(GetHitsCollection(CHCID, event)); 

//        std::cout << "eDep Fiber: " << eDep << std::endl;

        std::cout << "Population Fiber: " << pop << std::endl;

        G4AnalysisManager* analysis = G4AnalysisManager::Instance();
        //analysis->FillH1(CHCID, pop);
        analysis->FillH2(i, eDep, pop);
    }***/


/***
//  PMT Hits
    for(int i=0; i<10; i++) {
       
        //std::cout << "name: " << name[i] << j << scorer << std::endl;
        scorer = "/Pop";
        std::string scorerName = "sipm_" + name[i] + scorer;

        int CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
        //std::cout << "Name: " << scorerName << "CHCID: " << CHCID << std::endl;
        auto pop = GetSum(GetHitsCollection(CHCID, event)); 

        scorer = "/eDep";
        scorerName = "pmt_" + name[i] + scorer;
        CHCID = G4SDManager::GetSDMpointer()->GetCollectionID(scorerName);
        //std::cout << "Name: " << scorerName << "CHCID: " << CHCID << std::endl;
        auto eDep = GetSum(GetHitsCollection(CHCID,event));
        std::cout << "-----------PMT " << i << std::endl;
        std::cout << "Num Photons: " << pop << std::endl;
        std::cout << "eDep: "  << eDep << std::endl;
        G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    //    analysis->FillH1(CHCID, pop);

    }

***/



  if (HCE) {


    } else {
        G4cout << "No HCE" << G4endl;
    }
  //print per event (modulo n)
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    //G4cout << "---> End of event: " << eventID << G4endl;     
   // PrintEventStatistics(absoEdep, absoTrackLength, gapEdep, scintTrackLength);
  }
}  

//....

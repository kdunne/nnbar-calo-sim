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
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(scintHitsCollectionID == -1) {
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
  }


}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();

    int CHCID1 = -1;
    int CHCID2 = -1;

    if (CHCID2<0) {
        CHCID2 = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");
    }
    if (CHCID1<0){
        CHCID1 = G4SDManager::GetSDMpointer()->GetCollectionID("pmtHitCollection");
    }
  
    NNbarHitsCollection* AbsHits   = 0;
    NNbarHitsCollection* pmtHits   = 0;

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
        G4double X0     = 0.;
        G4double Y0     = 0.;
        G4double Z0     = 0.;
        G4double Xpos   = 0.;
        G4double Ypos   = 0.;
        G4double Zpos   = 0.;
        G4double eDepAbs = 0.;
        G4int cerenkovCounter = 0;
        G4double ke_init = 0.;
        G4double gammaA_KE   = -99.;
        G4double gammaB_KE   = -99.;
        G4double angle = -99.;
	G4ThreeVector gammaA_pos  = G4ThreeVector(-99.,0,0);
	G4ThreeVector gammaA_vert = G4ThreeVector(-99.,0,0);
	G4ThreeVector gammaB_pos  = G4ThreeVector(-99.,0,0);
	G4ThreeVector gammaB_vert = G4ThreeVector(-99.,0,0);
	G4ThreeVector gammaA_dist = G4ThreeVector(-99.,0,0);
	G4ThreeVector gammaB_dist = G4ThreeVector(-99.,0,0);


	pmtHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID1));
	AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID2));

	if(AbsHits) {
 
	    hitCount = AbsHits->entries();

            G4double x_pos[hitCount] = {0.};
 
            G4cout << "Processing " << hitCount << " hits" << G4endl;
           
	    for (G4int h=0; h<hitCount; h++) {
	        ltime               = ((*AbsHits)[h]) -> GetLocalTime();
		parentID 	    = ((*AbsHits)[h]) -> GetParentID();
     	        proc                = ((*AbsHits)[h]) -> GetProcess();
		G4double time       = ((*AbsHits)[h]) -> GetTime(); 
	        G4String name       = ((*AbsHits)[h]) -> GetName();
	    	G4int trID          = ((*AbsHits)[h]) -> GetTrackID();
		G4int i             = -99;
	        G4double kinEn      = ((*AbsHits)[h]) -> GetKinEn();
	        G4double eDep       = ((*AbsHits)[h]) -> GetEdep();
                G4double vertex_KE  = ((*AbsHits)[h]) -> GetVertexKE();
                G4ThreeVector pos   = ((*AbsHits)[h]) -> GetPos();	
                G4ThreeVector vert  = ((*AbsHits)[h]) -> GetVert();
                G4double isLastStep = ((*AbsHits)[h]) -> GetIsLast();
               
                eDepAbs += eDep/CLHEP::MeV;

                //G4cout << "eDepAbs: " << eDepAbs << G4endl;
                //G4cout << "eDep: " << eDep << G4endl;
                //analysis->FillH2(3, sqrt(pow(Xpos,2) + pow(Ypos,2))/CLHEP::cm, eDep/CLHEP::MeV);

                //analysis->FillNtupleIColumn(0,trID);
		//analysis->FillNtupleIColumn(0, event->GetEventID());
                //analysis->FillNtupleDColumn(1,sqrt(pow(Xpos,2) + pow(Ypos,2))/CLHEP::cm);
                //analysis->FillNtupleDColumn(2,eDep); 
                //analysis->AddNtupleRow(); 
                
                if (trID == 1){
	            ke_init = vertex_KE;
                    //std::cout << "KE init: " << ke_init << std::endl;
                    //if (isLastStep) {
                        //G4cout << "Filling...." << G4endl;
                        //analysis->FillH1(3, trackl/CLHEP::cm);
		        //analysis->FillH1(12, time/CLHEP::ns);
                    //}
	        }
     
                if (trID == 1 && (isLastStep || kinEn == 0) ) {
                    
                    G4ThreeVector distance = G4ThreeVector(vert.getX() - pos.getX(), 
                                                           vert.getY() - pos.getY(), 
                                                           vert.getZ() - pos.getZ());

                    trackl = sqrt ( pow(distance.getX(),2) + pow(distance.getY(),2) + pow(distance.getZ(),2) );
                    //std::cout << "trackl: " << trackl << std::endl;
                }

                //if (parentID == 1 && name=="gamma" ){
                    //G4cout << "trackID: " << trID <<  " particle: " << name << "KE: " << vertex_KE << G4endl;
               // }

                if (trID==2 && gammaA_KE == -99){
                    gammaA_KE = vertex_KE;
                    //gammaA_pos = pos;
                    //gammaA_vert = vert;
                    //std::cout << "gammaA Pos: " << pos.getX()/CLHEP::cm << std::endl; 
                }
                
                if (trID==3 && gammaB_KE ==-99){
                    gammaB_KE = vertex_KE;
                    //gammaB_pos = pos;
                    //gammaB_vert = vert;
		 } 

                if (trID == 2 && (isLastStep || kinEn==0)){
                    gammaA_pos = pos;
                    gammaA_vert = vert;
                    gammaA_dist = G4ThreeVector(gammaA_vert.getX() - gammaA_pos.getX(), 
                                                              gammaA_vert.getY() - gammaA_pos.getY(), 
                                                              gammaA_vert.getZ() - gammaA_pos.getZ());

                    std::cout << "Distance gamma 1 Last Step Z: " << gammaA_dist.getZ()/ CLHEP::cm << std::endl;
                }

                if (trID == 3 && (isLastStep || kinEn==0)){
                    gammaB_pos = pos;
                    gammaB_vert = vert;
                    gammaB_dist = G4ThreeVector(gammaB_vert.getX() - gammaB_pos.getX(), 
                                                              gammaB_vert.getY() - gammaB_pos.getY(), 
                                                              gammaB_vert.getZ() - gammaB_pos.getZ());

                    std::cout << "Distance gamma 2 Last Step Z: " << gammaB_dist / CLHEP::cm << std::endl;

                }

          
                if (ltime == 0 && name == "opticalphoton" && proc == "Cerenkov"){
		    cerenkovCounter++;
                }
	   
                if (name != "opticalphoton") {
                    //analysis->FillH1(5, Xpos/CLHEP::cm);
                    //analysis->FillH1(6, Ypos/CLHEP::cm);
                }

                if (eDep >0 ){
                    //analysis->FillH3(0, Xpos, Ypos, eDep);
                } 



	    }

            //G4cout << "eDepAbs: " << eDepAbs << G4endl;
            //G4cout << "ke_init: " << ke_init << G4endl;
            //G4cout << "frac: " << eDepAbs/ ke_init << G4endl;

            if (angle == -99 && (gammaA_KE > -99 || gammaB_KE > -99)) {

                G4double dotProd = (gammaA_dist.getX() * gammaB_dist.getX()) 
				+ (gammaA_dist.getY() * gammaB_dist.getY())
				+ (gammaA_dist.getZ() * gammaB_dist.getZ());


                G4double normA = sqrt(pow(gammaA_dist.getX(),2) + pow(gammaA_dist.getY(),2) + pow(gammaA_dist.getZ(),2));
                G4double normB = sqrt(pow(gammaB_dist.getX(),2) + pow(gammaB_dist.getY(),2) + pow(gammaB_dist.getZ(),2));

                G4int ratio = dotProd/(normA*normB);
                //std::cout << "ratio: " << ratio << std::endl;

                if (ratio == -1){
                    angle = acos(-1) * 180./M_PI;;
                } else {
                    angle = abs(acos(dotProd / (normA*normB)) * 180./M_PI);
    
                }
               
/***
 if (angle<180.){
                std::cout << "dotProd: " << dotProd << std::endl;
                std::cout << "normA: " << normA << std::endl;
                std::cout << "normB: " << normB << std::endl;
                std::cout << "dotProd/(normA*normB): " << dotProd/(normA*normB) << std::endl;
                std::cout << "acos(): " << acos(dotProd/(normA*normB)) << std::endl; 
                std::cout << "Angle [degrees]: " << angle << std::endl;  

                }
***/
                analysis->FillH2(4,gammaA_KE/CLHEP::MeV, gammaB_KE/CLHEP::MeV);
                analysis->FillH1(11,gammaA_KE/CLHEP::MeV);
                analysis->FillH1(12,gammaB_KE/CLHEP::MeV);
                analysis->FillH1(13, angle);
            }


            analysis->FillH1(10, eDepAbs / ke_init);

            /***G4double moliere, sum, mean, var, std = 0.;
            for(i=0; i < hitCount; i++) { sum += x_pos[i];}
            mean = sum / hitCount;
            for(i=0; i < hitCount; i++) { var += pow(x_pos[i] - mean, 2);};
            std = sqrt(var/hitCount);***/
 
            //G4cout << "Mean: " << mean/CLHEP::cm << " cm" << G4endl;
            //G4cout << "Std: " << std/CLHEP::cm << " cm" << G4endl;

            //moliere = 1.65*std;
            //G4cout << "moliere: " << moliere/CLHEP::cm << " cm" << G4endl;
            //analysis->FillH1(8, moliere/CLHEP::cm);

            

	    
            if (eDepAbs>0){                  
                if(cerenkovCounter>0){
                    //G4cout << "Filling ceren histos with " << cerenkovCounter << " photons and " << eDepAbs/CLHEP::MeV << " MeV deposited." << G4endl;
                    analysis->FillH1(3, trackl/CLHEP::cm);
                    analysis->FillH2(1, cerenkovCounter, eDepAbs/CLHEP::MeV);  
                    analysis->FillH2(0, trackl/CLHEP::cm, cerenkovCounter);
                    analysis->FillH1(4, eDepAbs/CLHEP::MeV); 
                    analysis->FillH1(0, cerenkovCounter);     
                }
            }
            //G4cout << "Total Edep in lead-glass: " << eDepAbs/CLHEP::MeV << G4endl;
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

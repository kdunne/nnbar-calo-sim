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
#include "RunAction.hh"
#include "HistoManager.hh"
#include "detSD.hh"
#include "NNbarHit.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "Randomize.hh"
#include "CLHEP/Random/JamesRandom.h"
#include <iomanip>
//....

extern G4double event_number; 
extern G4ThreadLocal G4int local_event_number;

extern std::ofstream CV_outFile;

namespace {G4Mutex RunActionMutex = G4MUTEX_INITIALIZER;}

EventAction::EventAction(HistoManager *histo): 
    G4UserEventAction(),
	fHistoManager(histo),
    detHitsCollectionID(-1)
{
}

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
//.....

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
	if(detHitsCollectionID == -1) {
		detHitsCollectionID = pSDManager->GetCollectionID("detLV/detHitCollection");
	}
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
	
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	int detHC = -1;
	if (detHC<0) {detHC = G4SDManager::GetSDMpointer()->GetCollectionID("detLV/detHitCollection");}
	NNbarHitsCollection* detHits  = 0;
	
	if (HCE) {

		G4AutoLock lock(&RunActionMutex);

		fHistoManager->ClearEventVectors();
		detHits = (NNbarHitsCollection*)(HCE->GetHC(detHC));
		if (detHits) {
			G4int hitCount = detHits->entries();
			for (G4int h=0; h<hitCount; h++) {
				G4String name = ((*detHits)[h]) -> GetName();
				if (name != "opticalphoton"){
					G4double time     = ((*detHits)[h]) -> GetTime(); 
					G4int pid   = ((*detHits)[h]) -> GetPID(); 
					G4double ekin   = ((*detHits)[h]) -> GetKinEn(); 
					G4double xx = ((*detHits)[h]) -> GetPosX();
					G4double yy = ((*detHits)[h]) -> GetPosY();
					G4double zz = ((*detHits)[h]) -> GetPosZ();
					G4double pX = ((*detHits)[h]) -> GetPX();
					G4double pY = ((*detHits)[h]) -> GetPY();
					G4double pZ = ((*detHits)[h]) -> GetPZ();
			
					fHistoManager->FillDetVectors(pid, time, ekin, xx, yy, zz, pX, pY, pZ);
				}
			}
		}

	}

	else {G4cout << "No HCE" << G4endl;}
	auto eventID = event->GetEventID();
	auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {G4cout << "---> End of event: " << eventID << G4endl;}

}  

#include "EventAction.hh"
#include "Analysis.hh"
#include "HistoManager.hh"
#include "GenericSD.hh"
#include "NNbarHit.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "Randomize.hh"
#include <iomanip>
#include <math.h>

//....
namespace {G4Mutex RunActionMutex = G4MUTEX_INITIALIZER;}

EventAction::EventAction(HistoManager* histo): 
    G4UserEventAction(),fHistoManager(histo),
    scintHitsCollectionID(-1)
{}

//....

EventAction::~EventAction()
{}

//....

G4THitsMap<G4double>* EventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
  auto hitsCollection = static_cast<G4THitsMap<G4double>*>(event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
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
	     scintHitsCollectionID = pSDManager->GetCollectionID("ScintDet_Collection");
	}
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
	if(scintHitsCollectionID  < 0) {return;}
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();

	int S_CID=-1;

	G4String name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

	//  Scintillator Bar Hits
	if(S_CID<0) S_CID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintDet_Collection");
	NNbarHitsCollection* scintHits=0;

	if (HCE) {
		//G4AutoLock lock(&RunActionMutex);
		fHistoManager->ClearEventVectors();

		scintHits = (NNbarHitsCollection*)(HCE->GetHC(S_CID));
		if (scintHits) {
			G4int evno = event->GetEventID();
			G4int hitCount = scintHits->entries();
			for (G4int h=0; h<hitCount; h++) { 
				G4String pname = ((*scintHits)[h]) -> GetName();
				G4String dname = ((*scintHits)[h]) -> GetDetName();
				G4String proc = ((*scintHits)[h]) -> GetProcess();
				if (dname == "Scintillator"){
					if (pname != "opticalphoton"){
						G4int parentID = ((*scintHits)[h]) -> GetParentID();
						G4int trID = ((*scintHits)[h]) -> GetTrackID();
						G4int PID = ((*scintHits)[h]) -> GetParticleID();
						G4double x = ((*scintHits)[h]) -> GetPos().getX();
						G4double y = ((*scintHits)[h]) -> GetPos().getY();
						G4double z = ((*scintHits)[h]) -> GetPos().getZ();

						G4double vertx = ((*scintHits)[h]) -> GetVert().getX();
						G4double verty = ((*scintHits)[h]) -> GetVert().getY();
						G4double vertz = ((*scintHits)[h]) -> GetVert().getZ();

						G4double time = ((*scintHits)[h]) -> GetTime();
						G4int xid = ((*scintHits)[h]) -> GetXID();

						G4double kinEn = ((*scintHits)[h]) -> GetKinEn();
						G4double eDep = ((*scintHits)[h]) -> GetEdep();
						G4int nPhotons = ((*scintHits)[h])->GetPhotons();

						fHistoManager->FillScintVectors(evno, trID, PID, parentID, xid, time, kinEn, eDep, nPhotons, x, y, z, vertx, verty, vertz);
					}
				}
				else if (dname == "WLSFiber"){
					if (pname == "opticalphoton"){
						//G4cout << proc << G4endl;
						G4int trID = ((*scintHits)[h]) -> GetTrackID();
						G4int parentID = ((*scintHits)[h]) -> GetParentID();
						G4double x = ((*scintHits)[h]) -> GetPos().getX();
						G4double y = ((*scintHits)[h]) -> GetPos().getY();
						G4double z = ((*scintHits)[h]) -> GetPos().getZ();

						G4double vertx = ((*scintHits)[h]) -> GetVert().getX();
						G4double verty = ((*scintHits)[h]) -> GetVert().getY();
						G4double vertz = ((*scintHits)[h]) -> GetVert().getZ();

						G4double time = ((*scintHits)[h]) -> GetTime();
						G4int xid = ((*scintHits)[h]) -> GetXID();
						G4int procid=0;
						if(proc=="Scintillation") procid=1;
						if(proc=="OpWLS") procid=2;
						if(proc=="Cerenkow") procid=3;
						fHistoManager->FillFiberVectors(evno, trID, parentID, procid, xid, time, x, y, z, vertx, verty, vertz);
					}
				}
				else if (dname == "SiPM"){
					if (pname == "opticalphoton"){
						G4int trID = ((*scintHits)[h]) -> GetTrackID();
						G4int parentID = ((*scintHits)[h]) -> GetParentID();
						G4double time = ((*scintHits)[h]) -> GetTime();
						G4int xid = ((*scintHits)[h]) -> GetXID();
						fHistoManager->FillSiPMVectors(evno, trID, parentID, xid, time);
					}
				}
			}
		}
	}
	else {
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

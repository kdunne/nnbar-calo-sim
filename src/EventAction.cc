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
    genHitsCollectionID(-1),
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
	     genHitsCollectionID = pSDManager->GetCollectionID("GenDet_Collection");
	     scintHitsCollectionID = pSDManager->GetCollectionID("Scint_DetHitCollection");
	}
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
	if(scintHitsCollectionID  < 0) {return;}
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();

	int G_CID=-1;
	int S_CID=-1;

	G4String name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

	//  Scintillator Bar Hits
	if(G_CID<0) G_CID = G4SDManager::GetSDMpointer()->GetCollectionID("GenDet_Collection");
	if(S_CID<0) S_CID = G4SDManager::GetSDMpointer()->GetCollectionID("Scint_DetHitCollection");
	NNbarHitsCollection* genHits=0;
	NNbarHitsCollection* scintHits=0;

	if (HCE) {
		//G4AutoLock lock(&RunActionMutex);
		fHistoManager->ClearEventVectors();

		genHits = (NNbarHitsCollection*)(HCE->GetHC(G_CID));
		if (genHits) {
			G4int evno = event->GetEventID();
			G4int hitCount = genHits->entries();
			G4int checktrack=0;
			G4String prevname = "void";
			for (G4int h=0; h<hitCount; h++) { 
				G4String pname = ((*genHits)[h]) -> GetName();
				G4String dname = ((*genHits)[h]) -> GetDetName();
				G4String proc = ((*genHits)[h]) -> GetProcess();
				if (dname == "Scintillator"){
					if (pname != "opticalphoton"){
						G4int parentID = ((*genHits)[h]) -> GetParentID();
						G4int trID = ((*genHits)[h]) -> GetTrackID();
						G4int PID = ((*genHits)[h]) -> GetParticleID();
						G4double x = ((*genHits)[h]) -> GetPos().getX();
						G4double y = ((*genHits)[h]) -> GetPos().getY();
						G4double z = ((*genHits)[h]) -> GetPos().getZ();

						G4double vertx = ((*genHits)[h]) -> GetVert().getX();
						G4double verty = ((*genHits)[h]) -> GetVert().getY();
						G4double vertz = ((*genHits)[h]) -> GetVert().getZ();

						G4double time = ((*genHits)[h]) -> GetTime();
						G4int xid = ((*genHits)[h]) -> GetXID();

						G4double kinEn = ((*genHits)[h]) -> GetKinEn();
						G4double eDep = ((*genHits)[h]) -> GetEdep();
						G4int nPhotons = ((*genHits)[h])->GetPhotons();

						fHistoManager->FillScintVectors(evno, trID, PID, parentID, xid, time, kinEn, eDep, nPhotons, x, y, z, vertx, verty, vertz);
					}
					prevname=dname;
				}
				else if (dname == "WLSFiber"){
					if (pname == "opticalphoton"){
						//G4cout << proc << G4endl;
						G4int trID = ((*genHits)[h]) -> GetTrackID();
						G4int parentID = ((*genHits)[h]) -> GetParentID();
						G4double x = ((*genHits)[h]) -> GetPos().getX();
						G4double y = ((*genHits)[h]) -> GetPos().getY();
						G4double z = ((*genHits)[h]) -> GetPos().getZ();

						G4double vertx = ((*genHits)[h]) -> GetVert().getX();
						G4double verty = ((*genHits)[h]) -> GetVert().getY();
						G4double vertz = ((*genHits)[h]) -> GetVert().getZ();

						G4double time = ((*genHits)[h]) -> GetTime();
						G4int xid = ((*genHits)[h]) -> GetXID();
						G4int procid=0;
						if(proc=="Scintillation") procid=1;
						if(proc=="OpWLS") procid=2;
						if(proc=="Cerenkov") procid=3;
						if(checktrack!=trID && procid==2){
						//if(prevname!=dname){
							fHistoManager->FillFiberVectors(evno, trID, parentID, procid, xid, time, x, y, z, vertx, verty, vertz);
							checktrack=trID;
						}
					}
					prevname=dname;
				}
			}
		}
		scintHits = (NNbarHitsCollection*)(HCE->GetHC(S_CID));
		if (scintHits) {
			G4int evno = event->GetEventID();
			G4int hitCount = scintHits->entries();
			for (G4int h=0; h<hitCount; h++) { 
				G4String pname = ((*scintHits)[h]) -> GetName();
				if (pname == "opticalphoton"){
					G4int trID = ((*scintHits)[h]) -> GetTrackID();
					G4int parentID = ((*scintHits)[h]) -> GetParentID();
					G4double x = ((*scintHits)[h]) -> GetPos().getX();
					G4double y = ((*scintHits)[h]) -> GetPos().getY();
					G4double z = ((*scintHits)[h]) -> GetPos().getZ();

					G4double time = ((*scintHits)[h]) -> GetTime();
					G4double kinEn = ((*scintHits)[h]) -> GetKinEn();
					G4int xid = ((*scintHits)[h]) -> GetXID();
					fHistoManager->FillSiPMVectors(evno, trID, parentID, xid, time, kinEn, x, y, z);
				}
			}
		}
	}
	else {
		G4cout << "No HCE" << G4endl;
	}
	//print per event (modulo n)
	//auto eventID = event->GetEventID();
	//auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	//if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
	//	G4cout << "---> End of event: " << eventID << G4endl;     
	//	 PrintEventStatistics(absoEdep, absoTrackLength, gapEdep, scintTrackLength);
	//}
}  

//....

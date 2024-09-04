#include "EventAction.hh"
#include "Analysis.hh"
#include "HistoManager.hh"
#include "GenericSD.hh"
#include "NNbarHit.hh"
#include "DetectorConstruction.hh"

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
#include <cmath>

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
  auto hitsCollection = dynamic_cast<G4THitsMap<G4double>*>(event->GetHCofThisEvent()->GetHC(hcID));
  
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

G4bool IsRightScintillatorName(const G4String& name,const std::vector<G4String>& names)
{
	for(const auto& scintName : names)
		if(name=="Scintillator_" + scintName)
			return true;

	return false;
}


void EventAction::EndOfEventAction(const G4Event* event)
{  
	if(scintHitsCollectionID  < 0) {return;}
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();

	int G_CID=-1;
	int S_CID=-1;

	auto Scintillators = DetectorConstruction::GetPointer()->GetScintillatorNames();
	auto HolePlacement = DetectorConstruction::GetPointer()->GetHolesPlacement();

	//  Scintillator Bar Hits
	if(G_CID<0) G_CID = G4SDManager::GetSDMpointer()->GetCollectionID("GenDet_Collection");
	if(S_CID<0) S_CID = G4SDManager::GetSDMpointer()->GetCollectionID("Scint_DetHitCollection");
	NNbarHitsCollection* genHits=0;
	NNbarHitsCollection* scintHits=0;
	G4double FirstScintHitTime = 0.0;


	if (HCE) {
		//G4AutoLock lock(&RunActionMutex);
		fHistoManager->ClearEventVectors();

		genHits = (NNbarHitsCollection*)(HCE->GetHC(G_CID));
		if (genHits) {
			G4int evno = event->GetEventID();
			G4int hitCount = genHits->entries();
			G4int checktrack=0;
			G4int checktrackScint=0;
			G4int checktrackHole=0;
			G4String prevname = "void";
			G4bool phot0 = false;

			std::vector<std::vector<G4int>> scintHitsCount;
			std::vector<std::vector<G4int>> holeAHitsCount;
			std::vector<std::vector<G4int>> holeBHitsCount;
			std::vector<std::vector<G4int>> fiberAHitsCount;
			std::vector<std::vector<G4int>> fiberBHitsCount;

			std::vector<std::vector<G4double>> hitTimes;
			std::vector<std::vector<G4double>> hitPos;
			std::vector<std::vector<G4int>> hitPDG;

			// std::vector<std::vector<G4double>> hitkEn;
			// std::vector<std::vector<G4double>> hiteDep;


			for(int i = 0; i < Scintillators.size(); i++)
			{
				std::vector<G4int> vecI = {0};
				std::vector<G4int> vecIE = {};
				std::vector<G4double> vecT = {0.0};
				std::vector<G4double> vecD = {0.0,0.0,0.0};

				scintHitsCount.push_back(vecI);
				holeAHitsCount.push_back(vecI);
				holeBHitsCount.push_back(vecI);
				fiberAHitsCount.push_back(vecI);
				fiberBHitsCount.push_back(vecI);

				hitTimes.push_back(vecT);
				hitPos.push_back(vecD);
				hitPDG.push_back(vecIE);

				// hitkEn.push_back(vecD);
				// hiteDep.push_back(vecD);
			}

			for (G4int h=0; h<hitCount; h++) {
				auto hit = ((*genHits)[h]);
				G4String pname = hit -> GetName();
				G4String dname = hit -> GetDetName();
				G4String proc = hit -> GetProcess();
				G4String lvname = hit -> GetLVName();
				G4String solidname = hit -> GetSolidName();
				G4String nextlvname = hit -> GetNextLVName();
				G4String sdname = hit -> GetSDName();

				 //G4cout << "(genHits) " << hit -> GetTime() << "ns Particle: " << pname << " PID:" << event -> GetEventID() << " Detector: " << dname << " Process: " << proc << " LogicalVolume: " << lvname << " NextLogicalVolume: " << nextlvname << " Solid: " << solidname << "SDName: "<< sdname<<G4endl;
				if (IsRightScintillatorName(dname,Scintillators)) {
					if (pname != "opticalphoton") {
						G4int parentID = hit -> GetParentID();
						G4int trID = hit -> GetTrackID();
						G4int PID = hit -> GetParticleID();
						G4double x = hit-> GetPos().getX();
						G4double y = hit-> GetPos().getY();
						G4double z = hit -> GetPos().getZ();

						G4double vertx = hit -> GetVert().getX();
						G4double verty = hit -> GetVert().getY();
						G4double vertz = hit -> GetVert().getZ();

						G4double time = hit -> GetTime();
						G4int xid = hit -> GetXID();

						G4double kinEn = hit -> GetKinEn();
						G4double eDep = hit -> GetEdep();
						G4int nPhotons = hit -> GetPhotons();

//						G4cout << "trID: " << trID << " sdname: " << sdname << " photons: " << nPhotons << " particle: " << pname << " x: " << x << " y: " << y << " z: " << z << G4endl;


						for(int ii = 0; ii < Scintillators.size(); ii++)
						{

							if("ScintillatorLV_"+Scintillators[ii]==sdname)
							{

								if(sdname==nextlvname) {
									scintHitsCount[ii][0]+=nPhotons;
								}
								if(hitTimes[ii][0]==0.0 && hit->GetTime() > 1)
								{
									hitTimes[ii][0] = hit->GetTime();
									hitPos[ii][0] = x;
									hitPos[ii][1] = y;
									hitPos[ii][2] = z;

									hitPDG[ii].push_back(hit->GetParticleID());
									// hitkEn[ii].push_back(kinEn);
									// hiteDep[ii].push_back(eDep);
								}
							}
						}


						// G4cout << "PreVol: " << sdname << " PostVol: " << nextlvname << " DName: " << dname << G4endl;
						//fHistoManager->FillScintVectors(evno, trID, PID, parentID, xid, time, kinEn, eDep, nPhotons, x, y, z, vertx, verty, vertz);
					}
					prevname=dname;
				}
				else if (dname == "WLSFiber"){
					if (pname == "opticalphoton"){
						G4int trID = hit -> GetTrackID();
						G4int parentID = hit-> GetParentID();
						G4double x = hit -> GetPos().getX();
						G4double y = hit -> GetPos().getY();
						G4double z = hit -> GetPos().getZ();

						G4double vertx = hit -> GetVert().getX();
						G4double verty = hit -> GetVert().getY();
						G4double vertz = hit -> GetVert().getZ();

						G4double time = hit -> GetTime();
						G4int xid = hit -> GetXID();
						G4int partID = hit -> GetParticleID();

						G4int procid=0;
						if(proc=="Scintillation") procid=1;
						if(proc=="OpWLS") procid=2;
						if(proc=="Cerenkov") procid=3;

						if(checktrack!=trID && procid==2){
							// if(!phot0)
							// {

							for(int ii = 0; ii < Scintillators.size(); ii++)
							{
								if (sdname=="FiberALV_"+Scintillators[ii])
									fiberAHitsCount[ii][0]++;
								if (sdname=="FiberBLV_"+Scintillators[ii])
									fiberBHitsCount[ii][0]++;
							}

							// G4cout << "FillFiberVectors opticalphoton Hit #" << h << " TraceID: " << trID << " HitTime: " << time << G4endl;
							// fHistoManager->FillFiberVectors(evno, trID, parentID, procid, xid, time, x, y, z, vertx, verty, vertz);
							checktrack=trID;

								// phot0 = true;
							// }
						}
					}
					prevname=dname;
				}
				else if (dname == "Hole")
				{
					if (pname == "opticalphoton")
					{
						G4int trID = hit -> GetTrackID();
						if (checktrackHole!=hit->GetTrackID())
						{
							for(int ii = 0; ii < Scintillators.size(); ii++)
							{
								if (sdname==(HolePlacement==1 ? "HoleALV_Tube_" : "HoleALV_Box_")+Scintillators[ii])
									holeAHitsCount[ii][0]++;
								if (sdname==(HolePlacement==1 ? "HoleBLV_Tube_" : "HoleBLV_Box_")+Scintillators[ii])
									holeBHitsCount[ii][0]++;
							}
							checktrackHole=trID;
						}
					}
				}
				else if (dname == "Target")
				{
					// if(proc!="primary" && proc=="hadElastic")
					// 	G4cout << "Post Target Process: " << proc << " Event Num: " << event->GetEventID() << G4endl;


					auto process = 0;
					if(proc=="primary") process = 1;
					if(proc=="hadElastic") process = 2;
					if(proc=="protonInelastic") process = 3;
					if(proc=="hIoni") process = 4;
					if(proc=="compt") process = 5;

					fHistoManager->FillTargetVectors(process);

					// auto process = 0;
					// if(proc=="primary") process = 1;
					// if(proc=="hadElastic") process = 2;
					// if(proc=="protonInelastic") process = 3;
					// if(proc=="compt") process = 4;
     //
					// fHistoManager->FillTargetVectors(process);
					// G4cout << "Particle: " << pname << " Process: " << proc << " SensitiveDet: " << sdname << G4endl;
				}

				// G4cout << "EventID: " << event->GetEventID() << " dynamicScintname: " << dynamicScintname << G4endl;
			}

			for(int ii = 0; ii < Scintillators.size(); ii++) {
				fHistoManager->FillScintVectors(ii, scintHitsCount[ii][0], holeAHitsCount[ii][0], holeBHitsCount[ii][0],
				                                fiberAHitsCount[ii][0], fiberBHitsCount[ii][0], hitTimes[ii][0],
				                                hitPos[ii], hitPDG[ii]);
//			if (scintHitsCount[ii][0]) G4cout<<"EventID: " << event->GetEventID()<< " ii = "<<ii<<" scintHitsCount[ii][0] = "<<scintHitsCount[ii][0] <<G4endl;
			}
		}


		scintHits = (NNbarHitsCollection*)(HCE->GetHC(S_CID));
		if (scintHits) {
			G4int evno = event->GetEventID();
			G4int hitCount = scintHits->entries();


			std::vector<std::vector<G4int>> sipmA0phot;
			std::vector<std::vector<G4int>> sipmB0phot;
			std::vector<std::vector<G4int>> sipmA1phot;
			std::vector<std::vector<G4int>> sipmB1phot;
			std::vector<std::vector<G4double>> sipmA0t;
			std::vector<std::vector<G4double>> sipmB0t;
			std::vector<std::vector<G4double>> sipmA1t;
			std::vector<std::vector<G4double>> sipmB1t;

			for(int i = 0; i < Scintillators.size(); i++)
			{
				std::vector<G4int> vec = {0};
				sipmA0phot.push_back(vec);
				sipmB0phot.push_back(vec);
				sipmA1phot.push_back(vec);
				sipmB1phot.push_back(vec);

				std::vector<G4double> vecD = {};
				sipmA0t.push_back(vecD);
				sipmB0t.push_back(vecD);
				sipmA1t.push_back(vecD);
				sipmB1t.push_back(vecD);
			}

			for (G4int h=0; h<hitCount; h++) {
				G4String pname = ((*scintHits)[h]) -> GetName();
				G4String dname = ((*scintHits)[h]) -> GetSDName();
				G4String proc = ((*scintHits)[h]) -> GetProcess();
				G4String lvname = ((*scintHits)[h]) -> GetLVName();
				G4String solidname = ((*scintHits)[h]) -> GetSolidName();
				G4int trID = ((*scintHits)[h]) -> GetTrackID();
				G4String nextlvname = ((*scintHits)[h]) -> GetNextLVName();


				if (pname == "opticalphoton"){
					G4int parentID = ((*scintHits)[h]) -> GetParentID();
					G4double x = ((*scintHits)[h]) -> GetPos().getX();
					G4double y = ((*scintHits)[h]) -> GetPos().getY();
					G4double z = ((*scintHits)[h]) -> GetPos().getZ();

					G4double time = ((*scintHits)[h]) -> GetTime();
					G4double kinEn = ((*scintHits)[h]) -> GetKinEn();
					G4int xid = ((*scintHits)[h]) -> GetXID();
					G4int phot = ((*scintHits)[h]) -> GetPhotons();


					for(int ii = 0; ii < Scintillators.size(); ii++)
					{
						if(dname=="SiPMA0LV_"+Scintillators[ii])
						{
							sipmA0phot[ii][0]++;
							sipmA0t[ii].push_back(time);
						}
						if(dname=="SiPMA1LV_"+Scintillators[ii])
						{
							sipmA1phot[ii][0]++;
							sipmA1t[ii].push_back(time);
						}
						if(dname=="SiPMB0LV_"+Scintillators[ii])
						{
							sipmB0phot[ii][0]++;
							sipmB0t[ii].push_back(time);
						}
						if(dname=="SiPMB1LV_"+Scintillators[ii])
						{
							sipmB1phot[ii][0]++;
							sipmB1t[ii].push_back(time);
						}
					}
					// G4cout << "(scintHits) Particle: " << pname << " Detector: " << dname << " Process: " << proc << " LogicalVolume: " << lvname << " Solid: " << solidname << " TraceID: " << trID << G4endl;
				}
			}
			for(int ii = 0; ii < Scintillators.size(); ii++)
			{
				// G4cout << "photA0: " << sipmA0phot[ii][0] << " photA1: " << sipmA1phot[ii][0] << " photB0: " << sipmB0phot[ii][0] << " photB1: " << sipmB1phot[ii][0] << G4endl;
				fHistoManager->FillSiPMVectors(ii,sipmA0phot[ii][0],sipmA1phot[ii][0],sipmB0phot[ii][0],sipmB1phot[ii][0]);
				fHistoManager->FillSiPMVectors(ii,sipmA0t[ii],sipmA1t[ii],sipmB0t[ii],sipmB1t[ii]);
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

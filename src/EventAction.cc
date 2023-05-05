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
//#include "Analysis.hh"
#include "HistoManager.hh"
#include "detSD.hh"
#include "CVSD.hh"
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
    detHitsCollectionID(-1),
    //det2HitsCollectionID(-1),
    CVHitsCollectionID(-1)
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

//....

G4double EventAction::CalcEnergy(G4double edep, G4double w, G4double r)
{
	
	CLHEP::HepRandomEngine* engine = new CLHEP::HepJamesRandom();
	CLHEP::RandGaussQ rGauss(engine);
  	
	const G4double ppMeV=17400*0.56;
	const G4double parr[3] = {0.139, 110.68*mm, 0.0687};
	const G4double parw[2] = {0.233, 1802.9*mm};

	r = (r<10*cm) ? r : 10*cm;
	G4double np = floor(edep*ppMeV);
	np *= (parr[0]*exp(-r/parr[1])+parr[2]);
	np *= (parw[0]*exp(-w/parw[1]));
	np = rGauss.shoot(np,sqrt(np));
	np = (np<0) ? -np : np;

	return np/ppMeV;  
}  

G4double EventAction::CalcTime(G4double w, G4double v)
{
	
	CLHEP::HepRandomEngine* engine = new CLHEP::HepJamesRandom();
	CLHEP::RandGaussQ rGauss(engine);
  	
	const G4double c_scint=CLHEP::c_light/1.59; //from Kuraray data sheet
	G4double t=sqrt(w*w+v*v)/c_scint;
	t = rGauss.shoot(t,0.5*ns);
	t = (t<0) ? -t : t;
	return t;  
}  

//.....

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
	if(detHitsCollectionID == -1) {
		detHitsCollectionID = pSDManager->GetCollectionID("detectorLV/detHitCollection");
	}
//	if(det2HitsCollectionID == -1) {
//		det2HitsCollectionID = pSDManager->GetCollectionID("det2LV/detHitCollection");
//	}
	if(CVHitsCollectionID == -1) {
		CVHitsCollectionID = pSDManager->GetCollectionID("CVLV/CVHitCollection");
	}
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
	const G4int nPlanes=12;
	const G4int nBars=64;
	const G4double c_scint=CLHEP::c_light/1.59; //from Kuraray data sheet
	const G4double AttLen=3.0*m; //from Kuraray data sheet
	
	if(CVHitsCollectionID  < 0) {return;}
 
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	int detHC = -1;
	if (detHC<0) {detHC = G4SDManager::GetSDMpointer()->GetCollectionID("detectorLV/detHitCollection");}
	NNbarHitsCollection* detHits  = 0;
	
//	int det2HC = -1;
//	if (det2HC<0) {det2HC = G4SDManager::GetSDMpointer()->GetCollectionID("det2LV/detHitCollection");}
//	NNbarHitsCollection* det2Hits  = 0;

	int CVHC = -1;
	if (CVHC<0) {CVHC = G4SDManager::GetSDMpointer()->GetCollectionID("CVLV/CVHitCollection");}
	NNbarHitsCollection* CVHits  = 0;

	if (HCE) {

		G4AutoLock lock(&RunActionMutex);

		//std::cout<< "Thread ID " << G4Threading::G4GetThreadId() << ":: after run  " << local_event_number << std::endl;
		//if((local_event_number%100)==0){
		//	std::cout<< "Event " << local_event_number << std::endl;
		//}
	
		fHistoManager->ClearEventVectors();
		fHistoManager->ClearDigiVectors();
		detHits = (NNbarHitsCollection*)(HCE->GetHC(detHC));
		if (detHits) {
			hitCount = detHits->entries();
			for (G4int h=0; h<hitCount; h++) {
				name     = ((*detHits)[h]) -> GetName();
				if (name != "opticalphoton"){
					time     = ((*detHits)[h]) -> GetTime(); 
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

	//	det2Hits = (NNbarHitsCollection*)(HCE->GetHC(det2HC));
	//	if (det2Hits) {
	//		hitCount = det2Hits->entries();
	//		for (G4int h=0; h<hitCount; h++) {
	//			name     = ((*det2Hits)[h]) -> GetName();
	//			if (name != "opticalphoton"){
	//				time     = ((*det2Hits)[h]) -> GetTime(); 
	//				G4int pid   = ((*det2Hits)[h]) -> GetParentID(); 
	//				G4double ekin = ((*det2Hits)[h]) -> GetKinEn();
	//				G4double xx = ((*det2Hits)[h]) -> GetPosX();
	//				G4double yy = ((*det2Hits)[h]) -> GetPosY();
	//				G4double zz = ((*det2Hits)[h]) -> GetPosZ();
	//				G4double pX = ((*det2Hits)[h]) -> GetPX();
	//				G4double pY = ((*det2Hits)[h]) -> GetPY();
	//				G4double pZ = ((*det2Hits)[h]) -> GetPZ();

	//				fHistoManager->FillDet2Vectors(pid, time, ekin, xx, yy, zz, pX, pY, pZ);
	//			}
	//		}
	//	}


		CVHits = (NNbarHitsCollection*)(HCE->GetHC(CVHC));
		if (CVHits) {
			G4int cvpid[1000]={0};
			G4double dep[1000]={0};
			G4double u[1000]={0}, v[1000]={0}, w[1000]={0};
			hitCount = CVHits->entries();
			for (G4int h=0; h<hitCount; h++) {
				name     = ((*CVHits)[h]) -> GetName();
				if (name != "opticalphoton"){
					ltime    = ((*CVHits)[h]) -> GetLocalTime();
					//parentID = ((*CVHits)[h]) -> GetParentID();
					proc     = ((*CVHits)[h]) -> GetProcess();
					name     = ((*CVHits)[h]) -> GetName();
					time     = ((*CVHits)[h]) -> GetTime(); 
					trID     = ((*CVHits)[h]) -> GetTrackID();
					G4int cvbar = ((*CVHits)[h]) -> GetStave_ID();
					G4int plane = ((*CVHits)[h]) -> GetGroup_ID();
					G4int planedir = ((*CVHits)[h]) -> GetXID();
					G4int n = plane*nBars+cvbar;
					kinEn    = ((*CVHits)[h]) -> GetKinEn();
					eDep     = ((*CVHits)[h]) -> GetEdep();
					dep[n] += eDep;
					cvpid[n] = ((*CVHits)[h]) -> GetPID();
					trackl   = ((*CVHits)[h]) -> GetTrackLength();	
					G4double xx = ((*CVHits)[h]) -> GetPosX();
					G4double yy = ((*CVHits)[h]) -> GetPosY();
					G4double zz = ((*CVHits)[h]) -> GetPosZ();
					G4String vol_name = ((*CVHits)[h]) -> GetVolName();
					G4double pX = ((*CVHits)[h]) -> GetPX();
					G4double pY = ((*CVHits)[h]) -> GetPY();
					G4double pZ = ((*CVHits)[h]) -> GetPZ();
					if(plane==0){ 
						u[n]+=((yy-1610)*eDep); 
						v[n]+=((xx-1500+200*cvbar)*eDep); 
						w[n]+=(zz*eDep);
					}
					else if(plane==1){ 
						u[n]+=((yy-1630)*eDep); 
						v[n]+=((zz-1500+200*cvbar)*eDep); 
						w[n]+=(xx*eDep);
					}
					else if(plane==2){ 
						u[n]+=((yy+1610)*eDep); 
						v[n]+=((xx-1500+200*cvbar)*eDep); 
						w[n]+=(zz*eDep);
					}
					else if(plane==3){ 
						u[n]+=((yy+1630)*eDep); 
						v[n]+=((zz-1500+200*cvbar)*eDep); 
						w[n]+=(xx*eDep);
					}
					else if(plane==4){ 
						u[n]+=((xx-1610)*eDep); 
						v[n]+=((yy-1500+200*cvbar)*eDep); 
						w[n]+=(zz*eDep);
					}
					else if(plane==5){ 
						u[n]+=((xx-1630)*eDep); 
						v[n]+=((zz-1500+200*cvbar)*eDep); 
						w[n]+=(yy*eDep);
					}
					else if(plane==6){ 
						u[n]+=((xx+1610)*eDep); 
						v[n]+=((yy-1500+200*cvbar)*eDep); 
						w[n]+=(zz*eDep);
					}
					else if(plane==7){ 
						u[n]+=((xx+1630)*eDep); 
						v[n]+=((zz-1500+200*cvbar)*eDep); 
						w[n]+=(yy*eDep);
					}
					else if(plane==8){ 
						u[n]+=((zz-1610)*eDep); 
						v[n]+=((xx-1500+200*cvbar)*eDep); 
						w[n]+=(yy*eDep);
					}
					else if(plane==9){ 
						u[n]+=((zz-1630)*eDep); 
						v[n]+=((yy-1500+200*cvbar)*eDep); 
						w[n]+=(xx*eDep);
					}
					else if(plane==10){ 
						u[n]+=((zz+1610)*eDep); 
						v[n]+=((xx-1500+200*cvbar)*eDep); 
						w[n]+=(yy*eDep);
					}
					else if(plane==11){ 
						u[n]+=((zz+1630)*eDep); 
						v[n]+=((yy-1500+200*cvbar)*eDep); 
						w[n]+=(xx*eDep);
					}
					
					fHistoManager->FillCVVectors(local_event_number, trID, cvpid[n], bar, plane, planedir, 
							time, kinEn, eDep, trackl, xx, yy, zz, pX, pY, pZ);
				}
			}
			for (G4int ii=0;ii<nPlanes*nBars;ii++)
			{
				if (dep[ii]>1*keV)
				{
					u[ii] /= dep[ii];
					v[ii] /= dep[ii];
					w[ii] /= dep[ii];
					G4double w1=w[ii]+1.6*m;
					G4double w2=1.6*m-w[ii];
					G4double v1=v[ii]+5*cm;
					G4double v2=v[ii]-5*cm;
					G4double r1=sqrt(u[ii]*u[ii]+v1*v1);
					G4double r2=sqrt(u[ii]*u[ii]+v2*v2);
					
					G4double t1=CalcTime(w1,v1);
					G4double t2=CalcTime(w2,v1);
					G4double t3=CalcTime(w1,v2);
					G4double t4=CalcTime(w2,v2);
									
					G4double e1 = CalcEnergy(dep[ii], w1, r1);
					G4double e2 = CalcEnergy(dep[ii], w2, r1);
					G4double e3 = CalcEnergy(dep[ii], w1, r2);
					G4double e4 = CalcEnergy(dep[ii], w2, r2);
					
					G4double post = 0.5*(t1+t3-t2-t4)*c_scint;
					G4double pose = 0.5*log((e2+e4)/(e1+e3))*AttLen;

					if (e1>1.*keV && e2>1.*keV && e3>1.*keV && e4>1.*keV){
						fHistoManager->FillCVDigiVectors(cvpid[ii], ii%nBars, ii/nBars, dep[ii], u[ii], v[ii], w[ii],
								e1, e2, e3, e4, t1, t2, t3, t4, post, pose); 
					}
				} 
			}
		}
	}

	else {G4cout << "No HCE" << G4endl;}
	auto eventID = event->GetEventID();
	auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {G4cout << "---> End of event: " << eventID << G4endl;}

}  

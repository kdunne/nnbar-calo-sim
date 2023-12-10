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
#include "BarHit.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4HCofThisEvent.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4GenericMessenger.hh"

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
    samplingHitsCollectionID(-1),
    CVHitsCollectionID(-1),
	fWriteSampling(true),
	fWriteCV(true)
{
	DefineCommands();
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

//G4double EventAction::CalcEnergy(G4double edep, G4double w, G4double r)
//{
//	
//	CLHEP::HepRandomEngine* engine = new CLHEP::HepJamesRandom();
//	CLHEP::RandGaussQ rGauss(engine);
//  	
//	const G4double ppMeV=17400*0.56;
//	const G4double parr[3] = {0.139, 110.68*mm, 0.0687};
//	const G4double parw[2] = {0.2395, 1860.*mm};
//
//	r = (r<10*cm) ? r : 10*cm;
//	G4double np = floor(edep*ppMeV);
//	np *= (parr[0]*exp(-r/parr[1])+parr[2]);
//	np *= (parw[0]*exp(-w/parw[1]));
//	np = rGauss.shoot(np,sqrt(np));
//	np = (np<0) ? -np : np;
//
//	return np/ppMeV;  
//}  
//
//G4double EventAction::CalcTime(G4double w, G4double v)
//{
//	
//	CLHEP::HepRandomEngine* engine = new CLHEP::HepJamesRandom();
//	CLHEP::RandGaussQ rGauss(engine);
//  	
//	const G4double c_scint=CLHEP::c_light/1.59; //from Kuraray data sheet
//	G4double t=sqrt(w*w+v*v)/c_scint;
//	t = rGauss.shoot(t,0.5*ns);
//	t = (t<0) ? -t : t;
//	return t;  
//}  

//.....

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
	if(detHitsCollectionID == -1) {
		detHitsCollectionID = pSDManager->GetCollectionID("detectorLV/detHitCollection");
	}
	if(samplingHitsCollectionID == -1) {
		samplingHitsCollectionID = pSDManager->GetCollectionID("samplingLV/samplingHitCollection");
	}
	if(CVHitsCollectionID == -1) {
		CVHitsCollectionID = pSDManager->GetCollectionID("CVLV/CVHitCollection");
	}
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
	
	if(CVHitsCollectionID  < 0) {return;}
 
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	int detHC = -1;
	if (detHC<0) {detHC = G4SDManager::GetSDMpointer()->GetCollectionID("detectorLV/detHitCollection");}
	NNbarHitsCollection* detHits  = 0;
	
	int samplingHC = -1;
	if (samplingHC<0) {samplingHC = G4SDManager::GetSDMpointer()->GetCollectionID("samplingLV/samplingHitCollection");}
	NNbarHitsCollection* samplingHits  = 0;

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
			G4int hitCount = detHits->entries();
			for (G4int h=0; h<hitCount; h++) {
				G4String name     = ((*detHits)[h]) -> GetName();
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

		if(fWriteSampling){
			samplingHits = (NNbarHitsCollection*)(HCE->GetHC(samplingHC));
			if (samplingHits) {
				G4int hitCount = samplingHits->entries();
				for (G4int h=0; h<hitCount; h++) {
					G4String name     = ((*samplingHits)[h]) -> GetName();
					if (name != "opticalphoton"){
						G4double time     = ((*samplingHits)[h]) -> GetTime(); 
						G4int trID   = ((*samplingHits)[h]) -> GetTrackID(); 
						G4int pid   = ((*samplingHits)[h]) -> GetPID(); 
						G4double ekin = ((*samplingHits)[h]) -> GetKinEn();
						G4double xx = ((*samplingHits)[h]) -> GetPosX();
						G4double yy = ((*samplingHits)[h]) -> GetPosY();
						G4double zz = ((*samplingHits)[h]) -> GetPosZ();
						G4double pX = ((*samplingHits)[h]) -> GetPX();
						G4double pY = ((*samplingHits)[h]) -> GetPY();
						G4double pZ = ((*samplingHits)[h]) -> GetPZ();

						fHistoManager->FillSamplingVectors(trID, pid, time, ekin, xx, yy, zz, pX, pY, pZ);
					}
				}
			}
		}


		CVHits = (NNbarHitsCollection*)(HCE->GetHC(CVHC));
		if (CVHits) {
			G4int cvpid[nBars]={0};
			BarHit bar[nBars];
			for(int ii=0; ii<nPlanes; ii++){
				for(int jj=0; jj<nBarsPerPlane; jj++){
					G4int n = ii*nBarsPerPlane+jj;
					bar[n].SetGeom(ii,jj,bar_length/2.,bar_length,bar_width,bar_thickness);
				}
			}
			//G4double dep[nBars]={0};
			//G4double u[nBars]={0}, v[nBars]={0}, w[nBars]={0};
			G4int hitCount = CVHits->entries();
			for (G4int h=0; h<hitCount; h++) {
				G4String name     = ((*CVHits)[h]) -> GetName();
				if (name != "opticalphoton"){
					G4int ltime    = ((*CVHits)[h]) -> GetLocalTime();
					//parentID = ((*CVHits)[h]) -> GetParentID();
					G4String proc     = ((*CVHits)[h]) -> GetProcess();
					G4String name     = ((*CVHits)[h]) -> GetName();
					G4double time     = ((*CVHits)[h]) -> GetTime(); 
					G4int trID     = ((*CVHits)[h]) -> GetTrackID();
					G4int cvbar = ((*CVHits)[h]) -> GetStave_ID();
					G4int plane = ((*CVHits)[h]) -> GetGroup_ID();
					G4int planedir = ((*CVHits)[h]) -> GetXID();
					G4int n = plane*nBarsPerPlane+cvbar;
					G4double kinEn    = ((*CVHits)[h]) -> GetKinEn();
					G4double eDep     = ((*CVHits)[h]) -> GetEdep();
					//dep[n] += eDep;
					cvpid[n] = ((*CVHits)[h]) -> GetPID();
					G4double trackl   = ((*CVHits)[h]) -> GetTrackLength();	
					G4double xx = ((*CVHits)[h]) -> GetPosX();
					G4double yy = ((*CVHits)[h]) -> GetPosY();
					G4double zz = ((*CVHits)[h]) -> GetPosZ();
					//G4String vol_name = ((*CVHits)[h]) -> GetVolName();
					G4double pX = ((*CVHits)[h]) -> GetPX();
					G4double pY = ((*CVHits)[h]) -> GetPY();
					G4double pZ = ((*CVHits)[h]) -> GetPZ();

					if(fWriteCV){
						fHistoManager->FillCVVectors(local_event_number, trID, cvpid[n], cvbar, plane, planedir, 
								time, kinEn, eDep, trackl, xx, yy, zz, pX, pY, pZ);
					}
					switch(planedir){
						case 0:
							bar[n].AddHit(yy,xx,zz,eDep,time);
							break;
						case 1:
							bar[n].AddHit(yy,zz,xx,eDep,time);
							break;
						case 2:
							bar[n].AddHit(xx,yy,zz,eDep,time);
							break;
						case 3:
							bar[n].AddHit(xx,zz,yy,eDep,time);
							break;
						case 4:
							bar[n].AddHit(zz,xx,yy,eDep,time);
							break;
						case 5:
							bar[n].AddHit(zz,yy,xx,eDep,time);
							break;
						default:
							break;
					}
				}
			}
			for (G4int ii=0;ii<nBars;ii++)
			{
				if (bar[ii].GetEDep()>1*keV)
				{
					bar[ii].AnalyzeHits();
					if (bar[ii].GetE1()>1.*keV && bar[ii].GetE2()>1.*keV && bar[ii].GetE3()>1.*keV && bar[ii].GetE4()>1.*keV){
						fHistoManager->FillCVDigiVectors(cvpid[ii], bar[ii].GetBar(), bar[ii].GetPlane(),
								bar[ii].GetEDep(), bar[ii].GetTime(), bar[ii].GetX(), bar[ii].GetY(), bar[ii].GetZ(),
								bar[ii].GetE1(), bar[ii].GetE2(), bar[ii].GetE3(), bar[ii].GetE4(),
								bar[ii].GetT1(), bar[ii].GetT2(), bar[ii].GetT3(), bar[ii].GetT4(), 
								bar[ii].GetPosT(), bar[ii].GetPosE()); 
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

void EventAction::DefineCommands()
{
	// Define /B5/generator command directory using generic messenger class
	fMessenger = new G4GenericMessenger(this, "/eventWrite/", "Controls writing events to file.");

	auto& writeSamplingCmd = fMessenger->DeclareProperty("writeSampling", fWriteSampling);
	G4String guidance = "Write data from sampling layer to file.\n";
	writeSamplingCmd.SetGuidance(guidance);
	writeSamplingCmd.SetParameterName("writeSampling", true);
	writeSamplingCmd.SetDefaultValue("true");

	auto& writeCVCmd = fMessenger->DeclareProperty("writeCV", fWriteCV);
	guidance = "Write raw CV data to file.\n";
	writeCVCmd.SetGuidance(guidance);
	writeCVCmd.SetParameterName("writeCV", true);
	writeCVCmd.SetDefaultValue("true");
	
}

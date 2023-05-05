#include "NNbarRun.hh"
#include "HistoManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
//....

NNbarRun::NNbarRun(HistoManager *histo):
	G4Run(), fHistoManager(histo)
{}

NNbarRun::~NNbarRun() {}

void NNbarRun::RecordEvent(const G4Event *evt)
{
	G4Run::RecordEvent(evt);
	fHistoManager->FillTree();
} 

void NNbarRun::Merge(const G4Run* aRun)
{
  const NNbarRun* localRun = static_cast<const NNbarRun*>(aRun);
  G4Run::Merge(aRun);
} 

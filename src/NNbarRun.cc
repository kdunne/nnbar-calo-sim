#include "NNbarRun.hh"
//....

NNbarRun::NNbarRun():
	G4Run()
{}

NNbarRun::~NNbarRun() {}

void NNbarRun::Merge(const G4Run* aRun)
{
  const NNbarRun* localRun = static_cast<const NNbarRun*>(aRun);
  G4Run::Merge(aRun);
} 

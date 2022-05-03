#include "G4Run.hh"

class NNbarRun : public G4Run
{
	public:
		NNbarRun();
		virtual ~NNbarRun();
		virtual void RecordEvent(const G4Event*);
		virtual void Merge(const G4Run*);
	private:
  // data members                   
  G4int  fAbsoEdepHCID;
  G4int  fGapEdepHCID;
  G4int  fAbsoTrackLengthHCID;
  G4int  fCerenkovHCID;
  G4int  fScintTrackLengthHCID;
  G4int  scintHitsCollectionID;
  G4int  absHitsCollectionID;
  G4int  tubeHitsCollectionID;
  G4int  TPCHitsCollectionID;
  G4int  PMTHitsCollectionID;
  G4int  ShieldHitsCollectionID;
  G4int  SiliconHitsCollectionID;
 
  G4double x=0;
  G4double y=0;
  G4double z=0;
  
};

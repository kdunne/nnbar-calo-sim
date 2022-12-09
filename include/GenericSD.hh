#ifndef GenericSD_h
#define GenericSD_h 1

#include "NNbarHit.hh"

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class GenericSD : public G4VSensitiveDetector
{
public:
    GenericSD(G4String name);
    ~GenericSD();
    
    void Initialize(G4HCofThisEvent*);
    
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
    void EndOfEvent(G4HCofThisEvent*HCE);
    
private:
    NNbarHitsCollection *HitsCollection;
    G4String sensitiveDetectorName;
};
#endif



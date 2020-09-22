#ifndef ScintillatorSD_h
#define ScintillatorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include "NNbarHit.hh"
//#include "ScintillatorHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ScintillatorSD : public G4VSensitiveDetector
{
public:
    ScintillatorSD(G4String name);
    ~ScintillatorSD();
    
    
    std::ofstream ofs;
    void Initialize(G4HCofThisEvent*);
    
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
    void EndOfEvent(G4HCofThisEvent*HCE);
    
private:
    NNbarHitsCollection *HitsCollection;
    G4String sensitiveDetectorName;
};
#endif



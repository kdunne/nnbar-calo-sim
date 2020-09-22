#ifndef StackingAction_H
#define StackingAction_H 1

#include "globals.hh"
#include "RunAction.hh"
#include "G4UserStackingAction.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern std::vector<std::vector<int>> scint_photon_collection;
extern std::vector<int> lead_glass_photon_all;
extern std::vector<int> lead_glass_photon_pri;

class StackingAction : public G4UserStackingAction
{
 
public:
	StackingAction();
	~StackingAction();

public:
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
	void NewStage();
	void PrepareNewEvent();
	
private:
	
	G4int cerenkovCounter_all;	G4int cerenkovCounter_pri;
	G4int scintCounter; // no longer useful 
	std::vector<int> scint_layer_photon{0,0,0,0,0,0,0,0,0,0}; // not used anymore! 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

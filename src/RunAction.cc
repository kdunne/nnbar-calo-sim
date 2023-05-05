#include "RunAction.hh"
//#include "Analysis.hh"
#include "HistoManager.hh"
#include "NNbarRun.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include <iostream>
#include <fstream>
#include <string>
//....

extern G4int run_number; 
extern G4double event_number;
extern std::vector<std::vector<G4double>> particle_gun_record;

RunAction::RunAction(HistoManager *histo)
 : G4UserRunAction(), fHistoManager(histo), fRun(0)
{ 
    fMessenger = new G4GenericMessenger(this, "/particle_generator/", "Name the particle for the file name");
}

//....

RunAction::~RunAction()
{}

//.....

G4Run* RunAction::GenerateRun()
{ 
	G4cout << "Generating run "<< G4endl; 
	return new NNbarRun(fHistoManager); 
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
	//inform the runManager to save random number seed
	//G4RunManager::GetRunManager()->SetRandomNumberStore(true);

	fHistoManager->Book();
	
}



void RunAction::EndOfRunAction(const G4Run* run)
{
	fHistoManager->Save();
	// write the PMT records 

	if(IsMaster())
	{
		// clear all stored data 
		event_number = 0.0;
		particle_gun_record.clear();
		run_number++;
	}

}

//....

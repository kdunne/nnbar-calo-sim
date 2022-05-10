#include "ActionInitialization.hh"
#include "config.h"

#include "PrimaryGeneratorAction.hh"
#include "G4MCPLGenerator.hh"
#include "PrimaryGeneratorAction_CRY.hh"   //CRY is not specified => we will use the normal particle gun :)

#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

//....

ActionInitialization::ActionInitialization()
	: G4VUserActionInitialization()
{}

//....

ActionInitialization::~ActionInitialization()
{;}

//....

void ActionInitialization::BuildForMaster() const
{
	SetUserAction(new RunAction);
}

//....

void ActionInitialization::Build() const
{ 

	SetUserAction(new PrimaryGeneratorAction);  

#if CRY_BUILD==1
	std::cout << "check .... CRY build is activated "<<std::endl;
	SetUserAction(new PrimaryGeneratorAction_CRY);
#endif

#if VERSION_MCPL ==1
	std::cout << " (((((( ****** ((()))))))  version MCPL" << std::endl;
	SetUserAction(new G4MCPLGenerator("./mcpl_files/NNBAR_mfro_signal_GBL_jbar_50k_9001.mcpl"));
#endif

	SetUserAction(new RunAction);
	EventAction* eventAction = new EventAction();
	SetUserAction(eventAction);
	SetUserAction(new SteppingAction());
}  

//....

#include "ActionInitialization.hh"

#include "HistoManager.hh"

#include "G4MCPLGenerator.hh"

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
	HistoManager *histo = new HistoManager();
	SetUserAction(new RunAction(histo));
}

//....

void ActionInitialization::Build() const
{ 

	HistoManager *histo = new HistoManager();
	SetUserAction(new G4MCPLGenerator(histo));

	RunAction* runAction = new RunAction(histo);
	SetUserAction(runAction);

	EventAction* eventAction = new EventAction(histo);
	SetUserAction(eventAction);
	SetUserAction(new SteppingAction());
}  

//....

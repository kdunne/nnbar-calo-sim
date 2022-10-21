#include "G4Run.hh"

class HistoManager;

class NNbarRun : public G4Run
{
	public:
		NNbarRun(HistoManager *histo);
		virtual ~NNbarRun();
		virtual void RecordEvent(const G4Event*);
		virtual void Merge(const G4Run*);
	private:
		HistoManager* fHistoManager;
  // data members                   
 
};

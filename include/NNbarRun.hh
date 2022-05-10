#include "G4Run.hh"

class NNbarRun : public G4Run
{
	public:
		NNbarRun();
		virtual ~NNbarRun();
		virtual void Merge(const G4Run*);
	private:
  // data members                   
 
};

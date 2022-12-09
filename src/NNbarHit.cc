#include "NNbarHit.hh"

//**********************MT
G4ThreadLocal G4Allocator<NNbarHit>* NNbarHitAllocator=0;
//**********************MT

NNbarHit::NNbarHit()
	: G4VHit()
{
	localTime = 0.;
	parentID = 0;
	process = "";
	time = 0.;
	name = "";
	particleID = 0;
	trackID = 0;
	isLast = false;
	xHitID = 0;
	pos = G4ThreeVector(0., 0., 0.);
	vert = G4ThreeVector(0., 0., 0.);
	energyDeposit = 0.;
	vertex_KE = 0.;
	kinEnergy = 0.;
	photons = 0;

}

NNbarHit::~NNbarHit()
{}

NNbarHit::NNbarHit(const NNbarHit& right)
	: G4VHit()
{
	localTime = right.localTime;
	parentID = right.parentID;
	process = right.process;
	time = right.time;
	name = right.name;
	particleID = right.particleID;
	trackID = right.trackID;
	isLast = right.isLast;
	xHitID = right.xHitID;
	pos = right.pos;
	vert = right.vert;
	energyDeposit = right.energyDeposit;
	vertex_KE = right.vertex_KE;
	kinEnergy = right.kinEnergy;
	photons = right.photons;
}

const NNbarHit& NNbarHit::operator=(const NNbarHit& right)
{
	localTime = right.localTime;
	parentID = right.parentID;
	process = right.process;
	time = right.time;
	name = right.name;
	particleID = right.particleID;
	trackID = right.trackID;
	isLast = right.isLast;
	xHitID = right.xHitID;
	pos = right.pos;
	vert = right.vert;
	energyDeposit = right.energyDeposit;
	vertex_KE = right.vertex_KE;
	kinEnergy = right.kinEnergy;
	photons = right.photons;

	return *this;
}

G4bool NNbarHit::operator==(const NNbarHit& right) const
{
	return(xHitID==right.xHitID);
}

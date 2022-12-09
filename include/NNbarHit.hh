#ifndef NNbarHit_h
#define NNbarHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class NNbarHit : public G4VHit
{
	public:
		NNbarHit();
		NNbarHit(const NNbarHit&);
		virtual ~NNbarHit();


		const NNbarHit& operator=(const NNbarHit&);

		G4bool operator==(const NNbarHit&) const;

		//******************************MT
		inline void* operator new(size_t);
		inline void operator delete(void*);
		//******************************MT

	private:
		G4double localTime;
		G4int parentID;
		G4String process;
		G4double time;
		G4String name;
		G4int particleID; 
		G4int trackID; 
		G4int xHitID; // Hit x voxel 
		G4bool isLast;
		G4ThreeVector pos;
		G4ThreeVector vert;
		G4double energyDeposit; // Energy deposit associated with the hit
		G4double vertex_KE;
		G4double kinEnergy;
		G4int photons;
		G4String detName;

	public:
		inline G4double GetLocalTime(){return localTime;}
		inline G4int GetParentID(){return parentID;}
		inline G4String GetProcess(){return process;}
		inline G4String GetName(){return name;}
		inline G4double GetTime(){return time;}
		inline G4int GetParticleID(){return particleID;}
		inline G4int GetTrackID(){return trackID;}
		inline G4int GetXID(){return xHitID;} // Get x index of the voxel 
		inline G4ThreeVector GetPos(){return pos;}
		inline G4ThreeVector GetVert(){return vert;}
		inline G4double GetEdep(){return energyDeposit;} // Get energy deposit
		inline G4double GetKinEn(){return kinEnergy;}
		inline G4double GetVertexKE(){return vertex_KE;}
		inline G4bool GetIsLast(){return isLast;}
		inline G4int GetPhotons(){return photons;}
		inline G4String GetDetName(){return detName;}

		inline void SetLocalTime(G4double ltime){localTime = ltime;}
		inline void SetParentID(G4int parent){parentID = parent;}
		inline void SetProcess(G4String p){process = p;}
		inline void SetName(G4String n){name = n;}
		inline void SetTime(G4double t){time = t;}
		inline void SetParticleID(G4int part){particleID = part;}
		inline void SetTrackID(G4int track){trackID = track;}
		inline void SetXID(G4int xID){xHitID = xID;}
		inline void SetPos(G4ThreeVector z){pos = z;}
		inline void SetVert(G4ThreeVector z){vert = z;}
		inline void SetEDep(G4double eDep){energyDeposit = eDep;}
		inline void SetKinEn(G4double kinEn){kinEnergy = kinEn;}
		inline void SetVertexKE(G4double vertKE){vertex_KE = vertKE;}
		inline void SetIsLast(G4bool last){isLast = last;}
		inline void SetPhotons(G4int ph){photons = ph;}
		inline void SetDetName(G4String dn){detName = dn;}

};

typedef G4THitsCollection<NNbarHit> NNbarHitsCollection;
//******************************MT
extern G4ThreadLocal G4Allocator<NNbarHit>* NNbarHitAllocator;
//******************************MT

inline void* NNbarHit::operator new(size_t){ 

	if(!NNbarHitAllocator) {
		NNbarHitAllocator= new G4Allocator<NNbarHit>;
	}

	void *aHit;

	aHit = (void *) NNbarHitAllocator->MallocSingle();
	return aHit;

}

inline void NNbarHit::operator delete(void *aHit) {
	if(!NNbarHitAllocator){ 
		NNbarHitAllocator= new G4Allocator<NNbarHit>;
	}

	NNbarHitAllocator->FreeSingle((NNbarHit*) aHit);
}

#endif

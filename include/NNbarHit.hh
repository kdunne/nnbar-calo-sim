//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef NNbarHit_h
#define NNbarHit_h 1


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"


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
    G4int trackID; 
    G4int xHitID; // Hit x voxel 
    G4double posZ;
    //G4int zHitID; // Hit z voxel
    //G4int yHitID; // Hit y voxel 
    G4double energyDeposit; // Energy deposit associated with the hit
    G4double kinEnergy;

public:
    // Methods to get the information - energy deposit and associated
    // position in the phantom - of the hits stored in the hits collection  
 
    inline G4int GetTrackID()
    {return trackID;}

    inline G4int GetXID() // Get x index of the voxel 
    {return xHitID;}

    //inline G4int GetZID() // Get y index of the voxel   
    //{return zHitID;}

    //inline G4int GetYID() // Get z index of the voxel  
    //{return yHitID;}

    inline G4double GetPosZ()
    {return posZ;}


    inline G4double GetEdep() // Get energy deposit
    {return energyDeposit;}
 
    inline G4double GetKinEn()
    {return kinEnergy;}

    // Methods to store the information of the hit ( energy deposit, position in the phantom )
    // in the hits collection

/***    inline void SetKinEnAndPosition(G4int xx, G4double zz, G4double eDep, G4double kinEn)
    {
        xHitID = xx;
        posZ = zz;
        //yHitID = yy;
        //zHitID = zz;
        energyDeposit = eDep;
        kinEnergy = kinEn;
  }***/

    inline void SetTrackID(G4int track){
        trackID = track;
    }


    inline void SetXID(G4int xID){
        xHitID = xID;
    }


     inline void SetPosZ(G4int z){
        posZ = z;
    }


    inline void SetEDep(G4int eDep){
        energyDeposit = eDep;
    }

    inline void SetKinEn(G4int kinEn){
        kinEnergy = kinEn;
    }


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

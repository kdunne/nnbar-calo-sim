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

#ifndef ScintillatorHit_h
#define ScintillatorHit_h 1


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"


class ScintillatorHit : public G4VHit
{
public:
  ScintillatorHit();
  ScintillatorHit(const ScintillatorHit&);
  virtual ~ScintillatorHit();
 
 
  const ScintillatorHit& operator=(const ScintillatorHit&);
 
  G4bool operator==(const ScintillatorHit&) const;

//******************************MT
inline void* operator new(size_t);
inline void operator delete(void*);
//******************************MT

private:
  G4int xHitID; // Hit x voxel 
  G4double posZ;
  //G4int zHitID; // Hit z voxel
  //G4int yHitID; // Hit y voxel 
  G4double energyDeposit; // Energy deposit associated with the hit
  G4double kinEnergy;

public:
  // Methods to get the information - energy deposit and associated
  // position in the phantom - of the hits stored in the hits collection  
 
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

  inline void SetKinEnAndPosition(G4int xx, G4double zz, G4double eDep, G4double kinEn)
  {
    xHitID = xx;
    posZ = zz;
    //yHitID = yy;
    //zHitID = zz;
    energyDeposit = eDep;
    kinEnergy = kinEn;
  }
};

typedef G4THitsCollection<ScintillatorHit> ScintillatorHitsCollection;
//******************************MT
extern G4ThreadLocal G4Allocator<ScintillatorHit>* ScintillatorHitAllocator;
//******************************MT

inline void* ScintillatorHit::operator new(size_t)
{
 
  
 if(!ScintillatorHitAllocator) 
  ScintillatorHitAllocator= new G4Allocator<ScintillatorHit>;
 void *aHit;

 aHit = (void *) ScintillatorHitAllocator->MallocSingle();
 return aHit;

}

inline void ScintillatorHit::operator delete(void *aHit)
{
if(!ScintillatorHitAllocator) 
  ScintillatorHitAllocator= new G4Allocator<ScintillatorHit>;

ScintillatorHitAllocator->FreeSingle((ScintillatorHit*) aHit);
}

#endif

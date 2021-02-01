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
//

#include "DetectorConstruction.hh"
#include "AbsorberSD.hh"
#include "pmtSD.hh"
#include "LightGuideSD.hh"
#include "WLSMaterials.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh" 
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4PSPopulation.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....
DetectorConstruction::~DetectorConstruction()
{ 
}

//....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....


void DetectorConstruction::DefineMaterials() { 

  //
  // -------------------------- Define Materials -------------------------------
  //


  // Use materials defined in NIST Manager LUT
  auto nistManager = G4NistManager::Instance();
  
  G4double a;
  G4double z;
  G4double density;

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

/***
  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elNa = nistManager->FindOrBuildElement("Na");
  G4Element* elCa = nistManager->FindOrBuildElement("Ca");
***/
  
  fMaterials = WLSMaterials::GetInstance();


}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {

  // 10  5x3x50 cm Scintillator Bar
  G4double WorldSizeX = 30.*cm;
  G4double WorldSizeY = 5.*cm;
  G4double WorldSizeZ = .5*m;
  G4double scintThickness = 3.*cm;

  G4double WLSfiberZ  = WorldSizeZ;
  G4double WLSfiberR  = 2.6*mm;

  G4double HoleRadius       = 2.9*mm;
  G4double HoleLength       = WLSfiberZ;
  G4double FiberRadius      = 0.5*mm;

  G4double WLSfiberOrigin = 0.0;

  // Get materials
  auto defaultMaterial    = G4Material::GetMaterial("Galactic");
 
  // World
  auto worldS 
    = new G4Box("World",           // its name
                 WorldSizeX/2., WorldSizeY/2., WorldSizeZ/2.); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  

  //--------------------------------------------------
  // Extrusion
  //--------------------------------------------------

  auto ExtrusionS =
        new G4Box("Extrusion", WorldSizeX/2, WorldSizeY/2, WorldSizeZ/2);

  auto ExtrusionLV =
        new G4LogicalVolume(ExtrusionS,
                            FindMaterial("Coating"),
                            "Extrusion");



  G4double fExtrusionReflectivity = 1.;
  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       fExtrusionReflectivity);

  G4MaterialPropertiesTable* TiO2SurfaceProperty =
                                             new G4MaterialPropertiesTable();

  G4double p_TiO2[] = {2.00*eV, 3.47*eV};
  const G4int nbins = sizeof(p_TiO2)/sizeof(G4double);

  G4double refl_TiO2[] = {fExtrusionReflectivity,fExtrusionReflectivity};
  assert(sizeof(refl_TiO2) == sizeof(p_TiO2));
  G4double effi_TiO2[] = {0, 0};
  assert(sizeof(effi_TiO2) == sizeof(p_TiO2));

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,nbins);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,nbins);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    ExtrusionLV,
                    "Extrusion",
                    worldLV,
                    false,
                    0);

  new G4LogicalSkinSurface("TiO2Surface",ExtrusionLV,TiO2Surface);

  std::string name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};
  G4double xPos[] = {-13.5, -10.5, -7.5, -4.5, -1.5, 1.5, 4.5, 7.5, 10.5, 13.5};
  G4double yPos[] = {0,0,0,0,0,0,0,0,0,0};
  G4double zPos[] = {0,0,0,0,0,0,0,0,0,0};

  for(int i=0; i<10; i++) {
//    for(int j=0; j<10; j++){
        std::cout << xPos[i] <<"," << yPos[i] << std::endl;
        //std::cout << "name: " << name[i] << j << std::endl;


        BuildScintBar(worldLV, 
                          xPos[i], 
                          yPos[i],
                          name[i],
                          zPos[i], 
                          scintThickness,
                          WorldSizeY,
                          WorldSizeZ);
   // }
  }

    
  return worldPV;
}

//....

void DetectorConstruction::BuildScintBar(G4LogicalVolume* logicExtrusion, G4double xPos, G4double yPos, std::string name, G4double zPos, G4double scintThickness, G4double WorldSizeY, G4double WorldSizeZ) {

  G4double HoleRadius   = 2.9*mm;
  G4double WLSfiberR    = 2.6*mm;
  G4double WLSfiberZ    = WorldSizeZ/2;


  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  auto scintDetector = new G4MultiFunctionalDetector(name);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);

  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------

  auto ScintillatorS 
    = new G4Box("Scintillator", scintThickness/2, WorldSizeY/2, WorldSizeZ/2);

  auto ScintillatorLV 
    = new G4LogicalVolume(
                 ScintillatorS,
                 FindMaterial("Polystyrene"),
                 name);

  auto ScintillatorPV 
    = new G4PVPlacement(
                 0,
                 G4ThreeVector(xPos*cm, yPos*cm, zPos*cm),
                 ScintillatorLV,
                 name,
                 logicExtrusion,
                 false,
                 0);

  //----------------------------------------------
  // Hole
  // ---------------------------------------------
  std::string hole_name = "hole" + name;

  auto HoleS 
    = new G4Tubs("Hole",
                 0.,
                 HoleRadius,
                 WorldSizeZ/2,
                 0.*deg,
                 360.*deg);

  auto HoleLV 
    = new G4LogicalVolume(
                 HoleS,
                 FindMaterial("G4_AIR"),
                 hole_name);

  auto holePV = new G4PVPlacement(0,
                    G4ThreeVector(0., 0., 0.),
                    HoleLV,
                    "Hole",
                    ScintillatorLV,
                    false,
                    0);

  //--------------------------------------------------
  // WLS Fiber
  //--------------------------------------------------
  std::string fiber_name = "fiber" + name;

  std::cout << "name: " << fiber_name << " Fiber size: " << WLSfiberR << " Fiber Xpos: " << xPos*cm << std::endl;



  auto FiberS 
    = new G4Tubs("WLSFiber",
                  0.,
                  WLSfiberR,
                  WorldSizeZ/2,
                  0.*deg,
                  360.*deg);

  auto FiberLV 
    = new G4LogicalVolume(FiberS,
                          FindMaterial("PMMA"),
                          fiber_name);

  //     logicWLSfiber->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
                                 
  auto FiberPV 
    = new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0.),
                        FiberLV,
                        "WLSFiber",
                        HoleLV,
                        false,
                        0);

}


void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
 

}


G4Material* DetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}


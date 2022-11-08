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

  fMaterials = WLSMaterials::GetInstance();

}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {

  // 10  5x3x50 cm scintillator bars
  G4double WorldSizeX = 30.*cm;

  // 5x3x50 cm scintillator bars   
  //G4double WorldSizeX = 40.*cm;

  G4double WorldSizeY = 5.*cm;
  G4double WorldSizeZ = 502*mm;

  G4double scintThickness = 3.*cm; 
  //G4double scintThickness = 4.*cm;

  G4double WLSfiberZ  = WorldSizeZ - 2*mm;
  G4double WLSfiberR  = 1.8*mm;

  G4double HoleRadius       = 1*mm;
  G4double HoleLength       = WLSfiberZ;
  G4double FiberRadius      = .7*mm;

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
  
// One bar at a time
// std::string name[] = {"A"};
//G4double xPos[] = {0.};
//G4double yPos[] = {0.};
//G4double zPos[] = {0.};


  // 10 Bars
  std::string name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

  // 3 cm thick bars
  G4double xPos[] = {-13.5, -10.5, -7.5, -4.5, -1.5, 1.5, 4.5, 7.5, 10.5, 13.5};

  // 4 cm thick bars
  //G4double xPos[] = {-18, -14, -10, -6, -2, 2, 6, 10, 14, 18};

  G4double yPos[] = {0,0,0,0,0,0,0,0,0,0};
  G4double zPos[] = {0,0,0,0,0,0,0,0,0,0};

  for(int i=0; i<10; i++) {
        //std::cout << xPos[i] <<"," << yPos[i] << std::endl;
        //std::cout << "name: " << name[i] << std::endl;

        BuildScintBar(worldLV, 
                          xPos[i], 
                          yPos[i],
                          name[i],
                          zPos[i], 
                          scintThickness,
                          WorldSizeY,
                          WorldSizeZ);
  }

    
  return worldPV;
}

//....

void DetectorConstruction::BuildScintBar(G4LogicalVolume* worldLV, G4double xPos, G4double yPos, std::string name, G4double zPos, G4double scintThickness, G4double WorldSizeY, G4double WorldSizeZ) {

  G4double HoleRadius   = 1*mm;
  G4double WLSfiberR    = .9*mm;
  G4double WLSfiberZ    = (WorldSizeZ-2*mm) /2.;
  //G4double WLSfiberZ  = WorldSizeZ/2;
  G4double CoatingThickness = .25*mm;
  G4double SiPMThickness = 2.*mm;

  G4double ScintZSize = WorldSizeZ-(SiPMThickness*2);
  G4double ScintXYSize = scintThickness-CoatingThickness;


  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);


  //--------------------------------------------------
  // Extrusion
  //--------------------------------------------------


  auto ExtrusionS =
        new G4Box("Extrusion", scintThickness/2, WorldSizeY/2, ScintZSize/2);

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
  G4int nbins = sizeof(p_TiO2)/sizeof(G4double);

  G4double refl_TiO2[] = {fExtrusionReflectivity,fExtrusionReflectivity};
  assert(sizeof(refl_TiO2) == sizeof(p_TiO2));
  G4double effi_TiO2[] = {0, 0};
  assert(sizeof(effi_TiO2) == sizeof(p_TiO2));

  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,nbins);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,nbins);

  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

  new G4PVPlacement(0,
                    G4ThreeVector(xPos*cm, yPos*cm, zPos*cm),
                    ExtrusionLV,
                    "Extrusion",
                    worldLV,
                    false,
                    0);

  new G4LogicalSkinSurface("TiO2Surface",ExtrusionLV,TiO2Surface);






  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------


 
  std::string scintName = "scint_" + name;

  auto ScintillatorS 
    = new G4Box("Scintillator", (scintThickness-CoatingThickness*2)/2, (WorldSizeY-CoatingThickness*2)/2, (ScintZSize-CoatingThickness*2)/2);

  auto ScintillatorLV 
    = new G4LogicalVolume(
                 ScintillatorS,
                 FindMaterial("Polystyrene"),
                 scintName);

  auto ScintillatorPV 
    = new G4PVPlacement(
                 0,
                 G4ThreeVector(0, 0, 0),
                 ScintillatorLV,
                 scintName,
                 ExtrusionLV,
                 false,
                 0);



 std::cout << "Built scintillator " << scintName << std::endl;


/***
  //----------------------------------------------
  // Hole 
  // ---------------------------------------------

  std::string hole_name = "hole" + name;

  auto HoleS 
    = new G4Tubs("Hole",
                 0.,
                 HoleRadius,
                 (ScintZSize-CoatingThickness*2)/2,
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

***/



//-----------------------------
// Hole A
// -----------------------------


  std::string hole_name = "holeA_" + name;

  auto HoleS
    = new G4Tubs("Hole",
                 0.,
                 HoleRadius,
                 (ScintZSize-CoatingThickness*2)/2,
                 0.*deg,
                 360.*deg);

  auto HoleLV_A
    = new G4LogicalVolume(
                 HoleS,
                 FindMaterial("G4_AIR"),
                 hole_name);

  auto HolePV_A = new G4PVPlacement(0,
                    G4ThreeVector(0., -1.25*cm, 0.),
                    HoleLV_A,
                    "Hole",
                    ScintillatorLV,
                    false,
                    0);


//-----------------------------
// Hole B
// -----------------------------


  hole_name = "holeB_" + name;

 /*** auto HoleS
    = new G4Tubs("Hole",
                 0.,
                 HoleRadius,
                 (ScintZSize-CoatingThickness*2)/2,
                 0.*deg,
                 360.*deg);
***/

  auto HoleLV_B
    = new G4LogicalVolume(
                 HoleS,
                 FindMaterial("G4_AIR"),
                 hole_name);

  auto HolePV_B = new G4PVPlacement(0,
                    G4ThreeVector(0., 1.25*cm, 0.),
                    HoleLV_B,
                    "Hole",
                    ScintillatorLV,
                    false,
                    0);







  //--------------------------------------------------
  // Cladding A
  //--------------------------------------------------

  auto cladS
    = new G4Tubs("Clad1",
                 WLSfiberR,
                 WLSfiberR + (0.042*mm/2),  // radius
                 (ScintZSize-CoatingThickness*2)/2,   // length
                 0.0*deg,
                 360*deg);

  auto cladLV_A
    = new G4LogicalVolume(cladS,
                          FindMaterial("Pethylene"),
                          "clad");

  auto cladPV_A 
    = new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0.),
                        cladLV_A,
                        "clad",
                        HoleLV_A,
                        false,
                        0);


  //--------------------------------------------------
  // Cladding B
  //--------------------------------------------------

  auto cladLV_B
    = new G4LogicalVolume(cladS,
                          FindMaterial("Pethylene"),
                          "clad");

  auto cladPV_B 
    = new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0.),
                        cladLV_B,
                        "clad",
                        HoleLV_B,
                        false,
                        0);




  //--------------------------------------------------
  // WLS Fiber A
  //--------------------------------------------------
  std::string fiberName = "fiberA_" + name;


  auto FiberS 
    = new G4Tubs("WLSFiber",
                  0.,
                  WLSfiberR,
                  (ScintZSize-CoatingThickness*2)/2,
                  0.*deg,
                  360.*deg);

  auto FiberLV_A 
    = new G4LogicalVolume(FiberS,
                          FindMaterial("PMMA"),
                          fiberName);

  //FiberLV->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
                                 
  auto FiberPV_A 
    = new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0.),
                        FiberLV_A,
                        "WLSFiber",
                        cladLV_A,
                        false,
                        0);



  //--------------------------------------------------
  // WLS Fiber B
  //--------------------------------------------------
  fiberName = "fiberB_" + name;


  auto FiberLV_B 
    = new G4LogicalVolume(FiberS,
                          FindMaterial("PMMA"),
                          fiberName);

  //FiberLV->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));
                                 
  auto FiberPV_B 
    = new G4PVPlacement(0,
                        G4ThreeVector(0.,0.,0.),
                        FiberLV_B,
                        "WLSFiber",
                        cladLV_B,
                        false,
                        0);









}




void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  std::string name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

  //auto scintName = "scint_A" ;


  for(int i=0; i<10; i++) {
    
 
    // declare scint as Primitive Scoreer
    auto scintName = "scint_" + name[i];

    auto scintDetector = new G4MultiFunctionalDetector(scintName);
    G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);

    G4VPrimitiveScorer* primitive = new G4PSPopulation("Pop");

    G4String fltName, particleName;
    G4SDParticleFilter* photonFilter =
        new G4SDParticleFilter(fltName="optPhoton", particleName="opticalphoton");
    primitive->SetFilter(photonFilter);
    scintDetector->RegisterPrimitive(primitive);

    SetSensitiveDetector(scintName, scintDetector);
    G4cout << "set sensitive detector " << scintName << " done." << G4endl;




    // declare fiberA as Primitive Scorer
    auto fiberName = "fiberA_" + name[i];

    auto fiberDetector = new G4MultiFunctionalDetector(fiberName);
    G4SDManager::GetSDMpointer()->AddNewDetector(fiberDetector);

    primitive = new G4PSPopulation("Pop");

    primitive->SetFilter(photonFilter);
    fiberDetector->RegisterPrimitive(primitive);

    SetSensitiveDetector(fiberName, fiberDetector);
    G4cout << "set sensitive detector " << fiberName << " done." << G4endl;


    // declare fiberB as Primitive Scorer
    fiberName = "fiberB_" + name[i];

    fiberDetector = new G4MultiFunctionalDetector(fiberName);
    G4SDManager::GetSDMpointer()->AddNewDetector(fiberDetector);

    primitive = new G4PSPopulation("Pop");

    primitive->SetFilter(photonFilter);
    fiberDetector->RegisterPrimitive(primitive);

    SetSensitiveDetector(fiberName, fiberDetector);
    G4cout << "set sensitive detector " << fiberName << " done." << G4endl;



  }

}


G4Material* DetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}


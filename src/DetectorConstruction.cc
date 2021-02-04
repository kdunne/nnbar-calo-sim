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

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSurface.hh"
#include "G4OpticalSurface.hh"

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


void DetectorConstruction::DefineMaterials()
{ 

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



  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elNa = nistManager->FindOrBuildElement("Na");
  G4Element* elCa = nistManager->FindOrBuildElement("Ca");


  // Lead-glass from PDG
  G4Material* Abs = new G4Material("Abs", 4.07*g/cm3, 5);
  Abs->AddElement(elO, 0.156453);
  Abs->AddElement(elSi, 0.080866);
  Abs->AddElement(elTi, 0.008092);
  Abs->AddElement(elAs, .002651);
  Abs->AddElement(elPb, 0.751938);

  // Optical Glass from G4
  G4Material* Glass = new G4Material("Glass", 2.4*g/cm3, 4);
  Glass->AddElement(elO, 0.4598  );
  Glass->AddElement(elNa, 0.096441);
  Glass->AddElement(elSi, 0.336553);
  Glass->AddElement(elCa, 0.107205);

  //
  // ----------------- Generate and Add Material Properties Table ----------------
  //


  /*********** Lead Glass ********************/
  G4MaterialPropertiesTable* absMPT = new G4MaterialPropertiesTable();

  // Datasheet gives photon wavelengths. Convert to Energies. Energies must be in Ascending order.
  G4double PhotonWavelength[] =
      { 2325.4*nm, 1970.1*nm, 1529.6*nm, 1060.0*nm,
        1014.0*nm, 852.10*nm, 706.50*nm, 656.30*nm,
        643.80*nm, 632.80*nm, 589.30*nm, 587.60*nm,
        546.10*nm, 486.10*nm, 480.00*nm, 435.80*nm,
        404.70*nm, 365.00*nm
       };

  const G4int nEntries = sizeof(PhotonWavelength)/sizeof(G4double);
  G4double PhotonEnergy[nEntries];
  for (int i=0; i < nEntries; ++i) {
    PhotonEnergy[i] = (1240.*nm/PhotonWavelength[i])*eV;
  };

  // Lead Glass Schott SF5
  G4double refractiveIndex[] =
        { 1.63289, 1.63785, 1.64359, 1.65104,
          1.65206, 1.65664, 1.66327, 1.66661,
          1.66756, 1.66846, 1.67252, 1.67270,
          1.67764, 1.68750, 1.68876, 1.69986,
          1.71069, 1.73056};

  absMPT->AddProperty("RINDEX", PhotonEnergy, refractiveIndex, nEntries)->SetSpline(true);
  Abs->SetMaterialPropertiesTable(absMPT);


  /*********** Light Guide ********************/
  G4MaterialPropertiesTable* lightGuideMPT = new G4MaterialPropertiesTable();

  // SCHOTT Datasheet same wavelengths for lead-glass.
  // Optical Glass Schott N-BK7
  G4double LGrefractiveIndex[] =
        { 1.48921, 1.49495, 1.50091, 1.50669,
          1.50731, 1.50980, 1.51289, 1.51432,
          1.51472, 1.51509, 1.51673, 1.51680,
          1.51872, 1.52238, 1.52283, 1.52668,
          1.53024, 1.53627};

  lightGuideMPT->AddProperty("RINDEX", PhotonEnergy, LGrefractiveIndex, nEntries)->SetSpline(true);
  Glass->SetMaterialPropertiesTable(lightGuideMPT);


  // Print materials
  //G4cout << "Absorber Properties -------" << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  //absMPT->DumpTable();
}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double absoThickness = 50.*cm;
  //G4double absoThickness = 35.*cm;
  //G4double absoThickness = 30.*cm;
  //G4double absoThickness = 25.*cm;
  //G4double absoThickness = 22.*cm;
  //G4double absoThickness = 20.*cm;
  //G4double absoThickness = 18.*cm;
  //G4double absoThickness = 15.*cm;
  //G4double absoThickness = 10.*cm;

  G4double rangeThickness = 30.*cm;
  G4double gap = 5.*cm;
  G4double layerThickness = 3.*cm;
  G4double lightGuideThickness = 10.*cm;
  G4double pmtThickness  = 4.*cm;
  G4double calorSizeXY   = 50.*cm;
  G4double worldSizeXY = 50.*cm;
  //G4double worldSizeXY = 25.*cm;
  //G4double worldSizeXY   = 1. * calorSizeXY;
  //G4double worldSizeZ    = rangeThickness + gap + absoThickness + lightGuideThickness + pmtThickness;
  G4double worldSizeZ    = rangeThickness + absoThickness + lightGuideThickness + pmtThickness;

  // Z Positions
  G4double absPos = 0.5 * (0.5*worldSizeZ - (absoThickness - worldSizeZ*0.5));
  //G4double pmtPos = (-0.5*worldSizeZ + 0.5*pmtThickness);
  
  G4double lightGuidePos = (absPos - absoThickness/2.) - (lightGuideThickness / 2.);
  G4double pmtPos = (lightGuidePos - lightGuideThickness/2.) - (pmtThickness/2.);

  // Get materials
  auto defaultMaterial    = G4Material::GetMaterial("Galactic");
  auto absorberMaterial   = G4Material::GetMaterial("Abs");
  auto lightGuideMaterial = G4Material::GetMaterial("Glass");

 
  if ( ! defaultMaterial || ! absorberMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  

   
  // World
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.); // its size
                         
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
  
  //BuildAbsorberGrid(worldLV, AbsorberSD)

  // Absorber

  G4double xGrid[] = {0};
  G4double yGrid[] = {0};

  //G4double xGrid[] = {-10., -5., 0., 5., 10.};
  //G4double yGrid[] = {-10., -5., 0., 5., 10.};

//  G4double xGrid[] = {-8., -4., 0., 4., 8.};
//  G4double yGrid[] = {-8., -4., 0., 4., 8.};

//  G4double xGrid[] = {-10., -8., -6., -4., -2., 2., 4., 6., 8., 10.};
//  G4double yGrid[] = {-10., -8., -6., -4., -2., 2., 4., 6., 8., 10.};


  G4String name[] = {"A", "B", "C", "D", "E"};
  //G4String name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};
  G4String LGname[] = {"lgA", "lgB", "lgC", "lgD", "lgE"};
  //G4String LGname[] = {"lgA", "lgB", "lgC", "lgD", "lgE", "lgF", "lgG", "lgH", "lgI", "lgJ"};
  G4String pmtname[] = {"pmtA", "pmtB", "pmtC", "pmtD", "pmtE"};
  //G4String pmtname[] = {"pmtA", "pmtB", "pmtC", "pmtD", "pmtE", "pmtF", "pmtG", "pmtH", "pmtI", "pmtJ"};

  std::cout << "----Grid Positions-----" << std::endl;
  for(int i=0; i<1; i++) {
    for(int j=0; j<1; j++){
        std::cout << xGrid[i] <<"," << yGrid[j] << std::endl;
//        std::cout << "name: " << name[i] << j << std::endl;

        std::string absoName = name[i] + std::to_string(j);
        std::string lgName = LGname[i] + std::to_string(j);
        std::string pmtName = pmtname[i] + std::to_string(j);
//        std::cout << "name: " << absoName << std::endl;

        BuildAbsorberGrid(worldLV, 
                          xGrid[i], 
                          yGrid[j],
                          absoName,
                          absPos, 
                          calorSizeXY, 
                          absoThickness);

        BuildLightGuideGrid(worldLV, 
                          xGrid[i], 
                          yGrid[j],
                          lgName,
                          lightGuidePos, 
                          calorSizeXY, 
                          lightGuideThickness);

        BuildPMTGrid(worldLV, 
                          xGrid[i], 
                          yGrid[j],
                          pmtName,
                          pmtPos, 
                          calorSizeXY, 
                          pmtThickness);



    }
  }

  






  // ----- Surfaces------
/***
  G4OpticalSurface* opAbsSurf = new G4OpticalSurface("AbsSurf");

  opAbsSurf->SetType(dielectric_dielectric);
  opAbsSurf->SetFinish(polished);
  opAbsSurf->SetModel(unified); 


  G4LogicalBorderSurface* AbsSurf =
    new G4LogicalBorderSurface("AbsSurf",    // its name
                                absorberPV,  // PhysicalVolume 1
                                worldPV,     // Physical Volume 2
                                opAbsSurf);  //surface Property

  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
        (AbsSurf->GetSurface(absorberPV,worldPV)->GetSurfaceProperty());
  if (opticalSurface) opticalSurface->DumpInfo();
***/




//  G4VisAttributes* pmtVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,2.0));
//  pmtVisAtt->SetVisibility(true);
//  pmtLV->SetVisAttributes(pmtVisAtt);

  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  worldVisAtt->SetVisibility(false);
  worldLV->SetVisAttributes(worldVisAtt);

 
  // print parameters
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << absoThickness/CLHEP::cm << " cm"
    << "------------------------------------------------------------" << G4endl;
 
//  G4cout << "Optical Surface Info:" << G4endl;
//  if (opticalSurface) { opticalSurface->DumpInfo(); }

    
  return worldPV;
}

//....



void DetectorConstruction::BuildAbsorberGrid(G4LogicalVolume* worldLV, G4double xPos, G4double yPos, std::string name, G4double zPos, G4double calorSizeXY, G4double absoThickness)
{

  // declare absorber as AbsorberSD 

  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //auto absorberDetector = new G4MultiFunctionalDetector(name);
  //G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);

  auto absorberMaterial   = G4Material::GetMaterial("Abs");

  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeXY/2., calorSizeXY/2., absoThickness/2.); // its size
                        
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 name);          // its name

  auto absorberPV  
    = new G4PVPlacement(
                 0,                // no rotation
		 G4ThreeVector(xPos*cm, yPos*cm, zPos),  // its position 
                 absorberLV,       // its logical volume                         
                 name,           // its name
                 worldLV,	   // its mother volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

//  std::cout << "pos: " << xPos << "," << yPos << " name: " << name << std::endl; 



  // declare absorber as AbsorberSD
//  AbsorberSD* absorberDetector = new AbsorberSD(name);
//  G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);
//  SetSensitiveDetector(name, absorberDetector);

  G4String absorberDetectorName = "abs_" + name ;
  auto absorberDetector = new G4MultiFunctionalDetector(absorberDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);


  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("eDep");
  absorberDetector->RegisterPrimitive(primitive);


  primitive = new G4PSPopulation("pop");
  G4String fltName, particleName;
  G4SDParticleFilter* photonFilter = 
	new G4SDParticleFilter(fltName="optPhoton", particleName="opticalphoton");
  primitive->SetFilter(photonFilter);
  absorberDetector->RegisterPrimitive(primitive);

  SetSensitiveDetector(name, absorberDetector);


  G4VisAttributes* absorberVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  absorberVisAtt->SetVisibility(true);
  absorberLV->SetVisAttributes(absorberVisAtt);

}

void DetectorConstruction::BuildLightGuideGrid(G4LogicalVolume* worldLV, G4double xPos, G4double yPos, std::string name, G4double zPos, G4double calorSizeXY, G4double lightGuideThickness) {

  auto lightGuideMaterial = G4Material::GetMaterial("Glass");


//  G4double dx1 = 3.*cm/2.;
//  G4double dx2 = 2.*cm/2.;
//  G4double dy1 = 3.*cm/2.;
//  G4double dy2 = 2.*cm/2.;
//  G4double dz = lightGuideThickness;

  // Light Guide

  G4double dx2 = 10.*cm/2.; 
  G4double dx1 = dx2/2.;
  G4double dy2 = 10*cm/2.;
  G4double dy1 = dy2/2.;
  G4double dz  = lightGuideThickness;


  auto lightGuideS 
    = new G4Box("lightGuide",
                 50*cm/2, // Dx1
                 50*cm/2,  // Dy2
                 5*cm/2); // Dz


 /*** auto lightGuideS 
    = new G4Trap("lightGuide",
                 dx1, // Dx1
                 dx2, // Dx2
                 dy1,      // Dy1
                 dy2,  // Dy2
                 dz/2); // Dz
***/
                         
  auto lightGuideLV
    = new G4LogicalVolume(
                 lightGuideS,        // its solid
                 lightGuideMaterial, // its material
                 "lightGuideLV");          // its name

  auto lightGuidePV  
    = new G4PVPlacement(
                 0,                // no rotation
		 G4ThreeVector(xPos*cm, yPos*cm, zPos),  // its position 
                 lightGuideLV,       // its logical volume                         
                 name,           // its name
                 worldLV,	   // its mother volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

}


void DetectorConstruction::BuildPMTGrid(G4LogicalVolume* worldLV, G4double xPos, G4double yPos, std::string name, G4double zPos, G4double calorSizeXY, G4double pmtThickness) {

  auto pmtMaterial = G4Material::GetMaterial("Glass");

  // PMT
  G4double innerRadius_pmt = 0.*cm;
  G4double outerRadius_pmt = 3.*cm/2.;
//  G4double outerRadius_pmt = calorSizeXY/4.;
  G4double height_pmt = pmtThickness;
  G4double startAngle_pmt = 0.*deg;
  G4double spanningAngle_pmt = 360.*deg;

  auto pmtS 
    = new G4Tubs("PMT",
                 innerRadius_pmt,
                 outerRadius_pmt,
                 height_pmt/2.,
                 startAngle_pmt,
                 spanningAngle_pmt);

                         
  auto pmtLV
    = new G4LogicalVolume(
                 pmtS,        // its solid
                 pmtMaterial, // its material
                 "pmtLV");          // its name

  auto pmtPV  
    = new G4PVPlacement(
                 0,                // no rotation
		 G4ThreeVector(xPos*cm, yPos*cm, zPos),  // its position 
                 pmtLV,       // its logical volume                         
                 name,           // its name
                 worldLV,	   // its mother volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

}

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
 

  //G4String absorberDetectorName = "AbsoLV" ;
  //AbsorberSD* absorberDetector = new AbsorberSD(absorberDetectorName);
  //G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);
  //SetSensitiveDetector("AbsoLV", absorberDetector);

  // declare light guide  as LightGuideSD
  G4String LightGuideDetectorName = "LightGuideLV" ;
  LightGuideSD* LightGuideDetector = new LightGuideSD(LightGuideDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(LightGuideDetector);
  //SetSensitiveDetector("lightGuideLV", LightGuideDetector);

  // declare PMT as Primitive Scoreer

  G4String pmtDetectorName = "pmtLV" ;
  auto pmtDetector = new G4MultiFunctionalDetector(pmtDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtDetector);

  G4VPrimitiveScorer* primitive;
  primitive = new G4PSPopulation("Pop");

  G4String fltName, particleName;
  G4SDParticleFilter* photonFilter = 
	new G4SDParticleFilter(fltName="optPhoton", particleName="opticalphoton");
  primitive->SetFilter(photonFilter);
  pmtDetector->RegisterPrimitive(primitive);

  //SetSensitiveDetector("pmtLV", pmtDetector);

}

//....

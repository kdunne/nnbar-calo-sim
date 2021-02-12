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
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "TubeSD.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

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
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  
  G4double a;
  G4double z;
  G4double density;


//--------------- Elements-----------------------------------

  // BC-408 taken from datasheet
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");

  G4Element* elAl = nistManager->FindOrBuildElement("Al"); 

  // Air
  G4Element* elN = nistManager->FindOrBuildElement("N");

  // Lead-glass
  // taken from PDG
  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  
  // TPC
  G4Element* elAr = nistManager->FindOrBuildElement("Ar");
 

//----------------- Materials -----------------------------

  // --------Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // --------Tube
  G4Material* Aluminum = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);

  // Silicon
  G4Material* Silicon = new G4Material("Silicon", z=14., a= 28.0855*g/mole, density = 2.33*g/cm3);

  // --------Air
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, 2); 
  Air->AddElement(elN, 0.7);
  Air->AddElement(elO, 0.3);

  // ---------FR4
  //----- Epoxy
  G4Material* Epoxy = new G4Material("Epoxy" , density=1.2*g/cm3, 2);
  Epoxy->AddElement(elH, 2);
  Epoxy->AddElement(elC, 2);
  //----- SiO2 (Quarz)
  G4Material* SiO2 = new G4Material("SiO2",density= 2.200*g/cm3, 2);
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO , 2);
  //FR4 (Glass + Epoxy)
  G4Material* FR4 = new G4Material("FR4" , density=1.86*g/cm3, 2);
  FR4->AddMaterial(Epoxy, 0.472);
  FR4->AddMaterial(SiO2, 0.528);
  //fr4Material = FR4;

  // ----------TPC
  // CO2
  G4Material* CO2 = new G4Material("CO2", density= 1.98*mg/cm3, 2);
  CO2->AddElement(elO, 2);
  CO2->AddElement(elC, 1);
  // Ar/CO2 80/20
  G4Material* Gas = new G4Material("Gas", density=1.3954*mg/cm3, 2);
  Gas->AddElement(elAr, .8);
  Gas->AddMaterial(CO2, .2);
  

  // ----------Scintillator
  G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
  Scint->AddElement(elH, 0.524573);
  Scint->AddElement(elC, 1 - 0.524573);

  // ----------Lead-glass Absorber
  G4Material* Abs = new G4Material("Abs", 3.86*g/cm3, 5);
  Abs->AddElement(elO, 0.156453);
  Abs->AddElement(elSi, 0.080866);
  Abs->AddElement(elTi, 0.008092);
  Abs->AddElement(elAs, .002651);
  Abs->AddElement(elPb, 0.751938);

//
// ----------------- Generate and Add Material Properties Table ----------------
//

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


//// Lead Glass Schott SF5

  G4double refractiveIndex[] =
        { 1.63289, 1.63785, 1.64359, 1.65104,
          1.65206, 1.65664, 1.66327, 1.66661,
          1.66756, 1.66846, 1.67252, 1.67270,
          1.67764, 1.68750, 1.68876, 1.69986,
          1.71069, 1.73056};


G4MaterialPropertiesTable* absMPT = new G4MaterialPropertiesTable();

absMPT->AddProperty("RINDEX", PhotonEnergy, refractiveIndex, nEntries)
      ->SetSpline(true);

Abs->SetMaterialPropertiesTable(absMPT);

// Print materials
  G4cout << "Absorber Properties -------" << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  absMPT->DumpTable();



//Scintillator Optical Properties

  /***
  const G4int nEntries2 = 12;

  G4double ScintPhotonEnergy[nEntries2] =

  { 2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV, 2.76*eV,
    2.82*eV, 2.92*eV, 2.95*eV, 3.02*eV, 3.1*eV,
    3.26*eV, 3.44*eV
  };

 //   {3.44*eV, 3.26*eV, 3.1*eV, 3.02*eV, 2.95*eV,
 //    2.92*eV, 2.82*eV, 2.76*eV, 2.7*eV, 2.58*eV,
 //    2.38*eV, 2.08*eV     
 //   };


  G4double rindex_scint[nEntries2] =
    {1.58, 1.58, 1.58, 1.58, 1.58,
     1.58, 1.58, 1.58, 1.58, 1.58,
     1.58, 1.58
    };

  G4double atten_scint[nEntries2] =
    {210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
     210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
     210*cm, 210*cm
    };

  G4double scintilFast[nEntries2] =
    {.04, .07, .20, .49, .84,
     1.0, .83, .55, .40, .17,
     .03, 0
    };

 G4double scintilSlow[nEntries2] =
    {.04, .07, .20, .49, .84,
     1.0, .83, .55, .40, .17,
     .03, 0
    };


  G4MaterialPropertiesTable *scintMPT = new G4MaterialPropertiesTable();
  scintMPT->AddProperty("RINDEX", ScintPhotonEnergy, rindex_scint, nEntries2)
	  ->SetSpline(true);
  scintMPT->AddProperty("ABSLENGTH", ScintPhotonEnergy, atten_scint, nEntries2)
	  ->SetSpline(true);
  scintMPT->AddProperty("FASTCOMPONENT", ScintPhotonEnergy, scintilFast, nEntries2)
	  ->SetSpline(true);
  scintMPT->AddProperty("SLOWCOMPONENT", ScintPhotonEnergy, scintilSlow, nEntries2)
	  ->SetSpline(true);

  // 64% of Antracene: 17400
  scintMPT->AddConstProperty("SCINTILLATIONYIELD", 11136000/keV);
  scintMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  scintMPT->AddConstProperty("FASTTIMECONSTANT", .9*ns);
  scintMPT->AddConstProperty("SLOWTIMECONSTANT", 2.1*ns);
  scintMPT->AddConstProperty("YIELDRATIO", 1.);

  Scint->SetMaterialPropertiesTable(scintMPT);

  // Print materials
  G4cout << "Scintillator Properties -------" << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  scintMPT->DumpTable();
***/
  
}
//....

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4int  nofLayers = 10;
  G4double  absoThickness  = 25.*cm;
  G4double  scintThickness =  3.*cm;
  G4double  calorSizeXY    =  98.*cm;
  G4double  tubeThickness  =  2.*cm;

  //auto calorThickness = (nofLayers*scintThickness) + absoThickness;
  auto calorThickness = (nofLayers*scintThickness) + absoThickness;
  //auto calorThickness = (nofLayers*scintThickness) + absoThickness + tubeThickness;
  auto worldSizeXY = 1 * calorSizeXY;
  auto worldSizeZ = 446*cm;
  //auto worldSizeZ  = 1 * calorThickness;


  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("Abs");
  auto scintMaterial = G4Material::GetMaterial("Scint");
  auto tubeMaterial = G4Material::GetMaterial("Aluminum"); 
  auto FR4Material = G4Material::GetMaterial("FR4");
  auto TPCMaterial = G4Material::GetMaterial("Gas");
  auto SiliconMaterial = G4Material::GetMaterial("Silicon");

 
  if ( ! defaultMaterial || ! absorberMaterial || ! scintMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  

   
  // World
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
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

  // Vacuum
  auto VacuumS 
    = new G4Box("Vacuum",           // its name
                 worldSizeXY/2, worldSizeXY/2, 112.06*cm/2); // its size
                         
  auto VacuumLV
    = new G4LogicalVolume(
                 VacuumS,           // its solid
                 defaultMaterial,  // its material
                 "Vacuum");         // its name
                                   
  auto VacuumPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 56.03*cm),  // at (0,0,0)
                 VacuumLV,          // its logical volume                         
                 "Vacuum",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  G4VisAttributes* VacuumVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  VacuumVisAtt->SetVisibility(true);
  VacuumLV->SetVisAttributes(VacuumVisAtt);


  // Silicon
  auto SiliconS 
    = new G4Box("Silicon",           // its name
                 worldSizeXY/2, worldSizeXY/2, .03*cm/2); // its size
                         
  auto FirstSiliconLV
    = new G4LogicalVolume(
                 SiliconS,           // its solid
                 SiliconMaterial,  // its material
                 "FirstSilicon");         // its name
                                   
  auto FirstSiliconPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 100.015*cm),  // at (0,0,0)
                 FirstSiliconLV,          // its logical volume                         
                 "FirstSilicon",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  auto SecondSiliconLV
    = new G4LogicalVolume(
                 SiliconS,           // its solid
                 SiliconMaterial,  // its material
                 "SecondSilicon");         // its name
                                   
  auto SecondSiliconPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 110.045*cm),  // at (0,0,0)
                 SecondSiliconLV,          // its logical volume                         
                 "SecondSilicon",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  G4VisAttributes* SiliconVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.));
  SiliconVisAtt->SetVisibility(true);
  FirstSiliconLV->SetVisAttributes(SiliconVisAtt);
  SecondSiliconLV->SetVisAttributes(SiliconVisAtt);



  /***
  // Calorimeter
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,    // its solid
                 defaultMaterial, // its material
                 "Calorimeter");  // its name
  
  auto calorPV  
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,0, 193.1*cm),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  *///
    
  // Absorber
  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");          // its name
 
  G4VisAttributes* absorberVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  absorberVisAtt->SetVisibility(true);
  absorberLV->SetVisAttributes(absorberVisAtt);
 

  auto absorberPV  
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 209.6*cm),
		 //G4ThreeVector(0., 0., 16.*cm),
                 //G4ThreeVector(0., 0., absoThickness/2. + (calorThickness/2. - absoThickness) ),
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 worldLV,	   // its mother volume
		 //layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  // Layer
  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, (scintThickness * nofLayers)/2); // its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 scintMaterial, //defaultMaterial,  // its material
                 "LayerLV");         // its name

  auto layerPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 181.*cm),
                 //G4ThreeVector(0., 0., -11.5*cm),
                 //G4ThreeVector(0., 0.,  -(scintThickness*nofLayers/2. ) + (calorThickness/2. - absoThickness) ), //  its position
                 layerLV,            // its logical volume                         
                 "Layer",            // its name
                 worldLV, // its mother volume
                 //calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

 /***
  // Scint
  auto scintS 
    = new G4Box("Scint",             // its name
                 calorSizeXY/2, calorSizeXY/2, scintThickness/2); // its size
                         
  auto scintLV
    = new G4LogicalVolume(
                 scintS,             // its solid
                 scintMaterial,      // its material
                 "ScintLV");      // its name
                                   
   auto scintPV
    = new G4PVReplica(
                 "Layer",          // its name
                 scintLV,          // its logical volume
                 layerLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 scintThickness);  // width of replica

***/

  // Vacuum Tube
  auto tubeS 
    = new G4Box("Tube",           // its name
                 calorSizeXY/2, calorSizeXY/2, tubeThickness/2); // its size
                         
  auto tubeLV
    = new G4LogicalVolume(
                 tubeS,           // its solid
                 tubeMaterial,  // its material
                 "TubeLV");         // its name

  auto tubePV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 113.06*cm),
                 //G4ThreeVector(0., 0.,  -27.5*cm), //  its position
                 tubeLV,            // its logical volume                         
                 "Tube",            // its name
                 worldLV, // its mother volume
                 //calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  G4VisAttributes* tubeVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  tubeVisAtt->SetVisibility(true);
  tubeLV->SetVisAttributes(tubeVisAtt);
 

  // Construct TPC
  //
  // FR4 front
  auto FR4frontS 
    = new G4Box("FR4front",           // its name
                 calorSizeXY/2, calorSizeXY/2, .16*cm/2.); // its size
                         
  auto FR4frontLV
    = new G4LogicalVolume(
                 FR4frontS,           // its solid
                 FR4Material,  // its material
                 "FR4frontLV");         // its name

  auto FR4frontPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 115.08*cm),
                 //G4ThreeVector(0., 0.,  -27.5*cm), //  its position
                 FR4frontLV,            // its logical volume                         
                 "FR4front",            // its name
                 worldLV, // its mother volume
                 //calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  // FR4 front
  auto FR4backS 
    = new G4Box("FR4bacl",           // its name
                 calorSizeXY/2, calorSizeXY/2, .16*cm/2.); // its size
                         
  auto FR4backLV
    = new G4LogicalVolume(
                 FR4backS,           // its solid
                 FR4Material,  // its material
                 "FR4backLV");         // its name

  auto FR4backPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 165.24*cm),
                 //G4ThreeVector(0., 0.,  -27.5*cm), //  its position
                 FR4backLV,            // its logical volume                         
                 "FR4back",            // its name
                 worldLV, // its mother volume
                 //calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 



  G4VisAttributes* FR4VisAtt= new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  FR4VisAtt->SetVisibility(true);
  FR4frontLV->SetVisAttributes(FR4VisAtt);
  FR4backLV->SetVisAttributes(FR4VisAtt); 


  // TPC
  auto TPCS 
    = new G4Box("TPC",           // its name
                 calorSizeXY/2, calorSizeXY/2, 50*cm/2.); // its size
                         
  auto TPCLV
    = new G4LogicalVolume(
                 TPCS,           // its solid
                 TPCMaterial,  // its material
                 "TPCLV");         // its name

  auto TPCPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 140.16*cm),
                 //G4ThreeVector(0., 0.,  -27.5*cm), //  its position
                 TPCLV,            // its logical volume                         
                 "TPC",            // its name
                 worldLV, // its mother volume
                 //calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  G4VisAttributes* TPCVisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  TPCVisAtt->SetVisibility(true);
  TPCLV->SetVisAttributes(TPCVisAtt);


  // print parameters
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << scintThickness/mm << "mm of " << scintMaterial->GetName() << " ] " << G4endl
    //<< tubeThickness/mm << "mm of " << tubeMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
 
    
  // Visualization attributes
  //worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  layerLV->SetVisAttributes(simpleBoxVisAtt);

  return worldPV;
}

//....

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare Scintillator as SinctillatorSD
  G4String scintDetectorName = "LayerLV" ;
  ScintillatorSD* scintDetector = new ScintillatorSD(scintDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);
  SetSensitiveDetector("LayerLV", scintDetector);
//  SetSensitiveDetector("ScintLV", scintDetector);

  // declare absorber as AbsorberSD
  G4String absorberDetectorName = "AbsoLV" ;
  AbsorberSD* absorberDetector = new AbsorberSD(absorberDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);
  SetSensitiveDetector("AbsoLV", absorberDetector);

  // declare vacuum as TubeSD
  G4String tubeDetectorName = "TubeLV" ;
  TubeSD* tubeDetector = new TubeSD(tubeDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(tubeDetector);
  SetSensitiveDetector("TubeLV", tubeDetector);


}

//....

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

G4ThreadLocal 
//G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//.....

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

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);


  // BC-408 taken from datasheet
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");

  G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
  Scint->AddElement(elH, 0.524573);
  Scint->AddElement(elC, 1 - 0.524573);

  // Lead-glass
  // taken from PDG

  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");


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
  G4double  absoThickness = 25.*cm;
  G4double  scintThickness =  3.*cm;
  G4double  calorSizeXY  = 1.*m;

  auto calorThickness = (nofLayers*scintThickness) + absoThickness;
  auto worldSizeXY = 1 * calorSizeXY;
  auto worldSizeZ  = 1 * calorThickness;


  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("Abs");
  auto scintMaterial = G4Material::GetMaterial("Scint");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! scintMaterial ) {
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
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
    
  // Absorber
  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");          // its name
 
 

  auto absorberPV  
    = new G4PVPlacement(
                 0,                // no rotation
		 G4ThreeVector(0., 0., absoThickness/2. + (calorThickness/2. - absoThickness) ),
                // G4ThreeVector(0., 0., -scintThickness/2), //  its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 calorLV,	   // its mother volume
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
                 defaultMaterial,  // its material
                 "LayerLV");         // its name

  auto layerPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,  -(scintThickness*nofLayers/2. ) + (calorThickness/2. - absoThickness) ), //  its position
                 layerLV,            // its logical volume                         
                 "Gap",            // its name
                 calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

 
  // Scint
  auto scintS 
    = new G4Box("Scint",             // its name
                 calorSizeXY/2, calorSizeXY/2, scintThickness/2); // its size
                         
  auto scintLV
    = new G4LogicalVolume(
                 scintS,             // its solid
                 scintMaterial,      // its material
                 "ScintLV");      // its name
                                   
 /***  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), //  its position
                 scintLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  ***/

   auto scintPV
    = new G4PVReplica(
                 "Layer",          // its name
                 scintLV,          // its logical volume
                 layerLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 scintThickness);  // width of replica


  // print parameters
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << scintThickness/mm << "mm of " << scintMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
 
    
  // Visualization attributes
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  // Scorers

  // declare Absorber as a MultiFunctionalDetector scorer
  auto absDetector = new G4MultiFunctionalDetector("Absorber");
  G4SDManager::GetSDMpointer()->AddNewDetector(absDetector);

  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("Edep");
  absDetector->RegisterPrimitive(primitive);


  // Hardcoded here for particle filter, must be changed for different particle
  primitive = new G4PSTrackLength("TrackLength");
  auto gun_part = new G4SDParticleFilter("gun_particle");
  gun_part->add("pi+");
  //gun_part->add("mu+");
  primitive->SetFilter(gun_part);
  absDetector->RegisterPrimitive(primitive);  

  primitive = new G4PSPopulation("Population");
  auto photonFilter = new G4SDParticleFilter("photonFilter");
  photonFilter->add("opticalphoton");
  primitive->SetFilter(photonFilter);
  absDetector->RegisterPrimitive(primitive);
  SetSensitiveDetector("AbsoLV",absDetector);
  
  // declare Scintillator as a MultiFunctionalDetector scorer
  auto scintDetector = new G4MultiFunctionalDetector("Scint");
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);

  primitive = new G4PSEnergyDeposit("Edep_gun");
  primitive->SetFilter(gun_part);
  scintDetector->RegisterPrimitive(primitive);
  
  primitive = new G4PSEnergyDeposit("Edep_other");
  primitive->SetFilter(gun_part);
  scintDetector->RegisterPrimitive(primitive);

  primitive = new G4PSTrackLength("TrackLength");
  primitive->SetFilter(gun_part);
  //primitive ->SetFilter(charged);
  scintDetector->RegisterPrimitive(primitive);  

  primitive = new G4PSPopulation("Population");
  primitive->SetFilter(photonFilter);
  scintDetector->RegisterPrimitive(primitive);

  SetSensitiveDetector("ScintLV", scintDetector);  

}

//....

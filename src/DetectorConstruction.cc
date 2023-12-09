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
#include "GenericSD.hh"
#include "Scint_DetSD.hh"
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
#include "G4GenericMessenger.hh"

//....

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),fLength(200.*cm),fWidth(20.*cm),fThickness(2.*cm)
{
	DefineCommands();
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

	G4double WorldSizeY = fWidth;
	G4double WorldSizeZ = fLength + 4*mm;

	//G4int scintBars = 10;
	G4double scintThickness = fThickness; 
	//G4double scintThickness = 4.*cm;

	G4double WLSfiberZ  = fLength;
	// G4double WLSfiberR  = 1.8*mm;

	// G4double HoleRadius       = 1*mm;
	// G4double HoleLength       = WLSfiberZ;
	// G4double FiberRadius      = .7*mm;

	// G4double WLSfiberOrigin = 0.0;

	// Get materials
	auto defaultMaterial    = G4Material::GetMaterial("Galactic");

	// World
	auto worldS = new G4Box("World", WorldSizeX/2., WorldSizeY/2., WorldSizeZ/2.);
	auto worldLV = new G4LogicalVolume(worldS, defaultMaterial, "World");
	auto worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "World", 0, false, 0, fCheckOverlaps);  

	// 10 Bars
	std::string name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

	//std::cout << "name: " << name[i] << std::endl;
	//G4double xPos = (i+0.5)*scintThickness-double(scintBars)/2*scintThickness;
	G4double xPos = 0.;
	G4double yPos = 0.;
	G4double zPos = 0.;

	BuildScintBar(worldLV, xPos, yPos, name[0], zPos, scintThickness, WorldSizeY, WorldSizeZ);
	return worldPV;
}

//....

void DetectorConstruction::BuildScintBar(G4LogicalVolume* worldLV, G4double xPos, G4double yPos, std::string name, G4double zPos, G4double scintThickness, G4double WorldSizeY, G4double WorldSizeZ) {

	G4double SiPMThickness = 2.*mm;
	G4double SiPMXYSize = 5.*mm;
	G4double HoleRadius   = 1*mm;
	G4double HoleY_A = -0.25*WorldSizeY;
	G4double HoleY_B = 0.25*WorldSizeY;
	G4double WLSfiberR    = 0.9*mm;
	G4double CoatingThickness = .25*mm;

	G4double ScintZSize = WorldSizeZ-(SiPMThickness*2);
	G4double ScintXYSize = scintThickness-CoatingThickness;

	G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

	//--------------------------------------------------
	// Extrusion
	//--------------------------------------------------

	auto ExtrusionS = new G4Box("Extrusion", scintThickness/2, WorldSizeY/2, ScintZSize/2);
	auto ExtrusionLV = new G4LogicalVolume(ExtrusionS, FindMaterial("Coating"), "Extrusion");

	G4cout << "Scintillator bar dimensions: " << fLength/cm << " cm x " << fWidth/cm << " cm x " << fThickness/cm << " cm" << G4endl;
	G4cout << "Checking G4Box size: " << ExtrusionS->GetZHalfLength() << " mm x " << ExtrusionS->GetYHalfLength() << " mm x " << ExtrusionS->GetXHalfLength() << " mm" << G4endl;
	
	G4double fExtrusionReflectivity = 1.;
	//G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface", glisur, ground, dielectric_metal, fExtrusionReflectivity);
	G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface");
	TiO2Surface->SetModel(unified);
	TiO2Surface->SetType(dielectric_metal);
	TiO2Surface->SetFinish(ground);

	G4MaterialPropertiesTable* TiO2SurfaceProperty = new G4MaterialPropertiesTable();

	std::vector<G4double> p_TiO2 = {2.00*eV, 3.47*eV};
	std::vector<G4double> refl_TiO2 = {fExtrusionReflectivity,fExtrusionReflectivity};
	std::vector<G4double> effi_TiO2 = {0, 0};
	TiO2SurfaceProperty->AddProperty("SPECULARLOBECONSTANT", p_TiO2, effi_TiO2);
	TiO2SurfaceProperty->AddProperty("SPECULARSPIKECONSTANT", p_TiO2, effi_TiO2);
	TiO2SurfaceProperty->AddProperty("BACKSCATTERCONSTANT", p_TiO2, effi_TiO2);
	TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2);
	TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2);

	TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

	new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), ExtrusionLV, "Extrusion", worldLV, false, 0);

	new G4LogicalSkinSurface("TiO2Surface",ExtrusionLV,TiO2Surface);

	TiO2Surface->DumpInfo();

	//--------------------------------------------------
	// Scintillator
	//--------------------------------------------------

	auto ScintillatorS = new G4Box("Scintillator", (scintThickness-CoatingThickness*2)/2, (WorldSizeY-CoatingThickness*2)/2, (ScintZSize-CoatingThickness*2)/2);

	auto ScintillatorLV = new G4LogicalVolume(ScintillatorS, FindMaterial("Polystyrene"), "ScintillatorLV");
	auto ScintillatorPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), ScintillatorLV, "Scintillator", ExtrusionLV, false, 0);

	//-----------------------------
	// Holes 
	// -----------------------------

	auto HoleS = new G4Tubs("Hole", 0., HoleRadius, (ScintZSize-CoatingThickness*2)/2, 0.*deg, 360.*deg);

	auto HoleLV_A = new G4LogicalVolume(HoleS, FindMaterial("G4_AIR"), "HoleALV");
	auto HolePV_A = new G4PVPlacement(0, G4ThreeVector(0., HoleY_A, 0.), HoleLV_A, "Hole", ScintillatorLV, false, 0);

	auto HoleLV_B = new G4LogicalVolume(HoleS, FindMaterial("G4_AIR"), "HoleBLV");
	auto HolePV_B = new G4PVPlacement(0, G4ThreeVector(0., HoleY_B, 0.), HoleLV_B, "Hole", ScintillatorLV, false, 1);

	//--------------------------------------------------
	// Cladding
	//--------------------------------------------------

	auto oCladS = new G4Tubs("oClad1", 0., WLSfiberR *1.06, (ScintZSize-CoatingThickness*2)/2, 0.0*deg, 360*deg);

	auto oCladLV_A = new G4LogicalVolume(oCladS, FindMaterial("FPethylene"), "oClad");
	auto oCladLV_B = new G4LogicalVolume(oCladS, FindMaterial("FPethylene"), "oClad");

	auto oCladPV_A = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), oCladLV_A, "oClad", HoleLV_A, false, 0);
	auto oCladPV_B = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), oCladLV_B, "oClad", HoleLV_B, false, 1);

	//   auto oCladPV_A = new G4PVPlacement(0, G4ThreeVector(0.,HoleY_A,0.), oCladLV_A, "oClad", ScintillatorLV, false, 0);
	//   auto oCladPV_B = new G4PVPlacement(0, G4ThreeVector(0.,HoleY_B,0.), oCladLV_B, "oClad", ScintillatorLV, false, 1);


	auto iCladS = new G4Tubs("iClad1", 0., WLSfiberR *1.03, (ScintZSize-CoatingThickness*2)/2, 0.0*deg, 360*deg);

	auto iCladLV_A = new G4LogicalVolume(iCladS, FindMaterial("Pethylene"), "iClad");
	auto iCladLV_B = new G4LogicalVolume(iCladS, FindMaterial("Pethylene"), "iClad");

	auto iCladPV_A = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), iCladLV_A, "iClad", oCladLV_A, false, 0);
	auto iCladPV_B = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), iCladLV_B, "iClad", oCladLV_B, false, 1);

	//  auto iCladPV_A = new G4PVPlacement(0, G4ThreeVector(0.,HoleY_A,0.), iCladLV_A, "iClad", ScintillatorLV, false, 0);
	//  auto iCladPV_B = new G4PVPlacement(0, G4ThreeVector(0.,HoleY_B,0.), iCladLV_B, "iClad", ScintillatorLV, false, 1);

	//--------------------------------------------------
	// WLS Fibers
	//--------------------------------------------------

	auto FiberS = new G4Tubs("WLSFiber", 0., WLSfiberR, (ScintZSize-CoatingThickness*2)/2, 0.*deg, 360.*deg);

	auto FiberLV_A = new G4LogicalVolume(FiberS, FindMaterial("PMMA"), "FiberALV");
	auto FiberLV_B = new G4LogicalVolume(FiberS, FindMaterial("PMMA"), "FiberBLV");

	//  auto FiberPV_A = new G4PVPlacement(0, G4ThreeVector(0.,HoleY_A,0.), FiberLV_A, "WLSFiber", ScintillatorLV, false, 0);
	//  auto FiberPV_B = new G4PVPlacement(0, G4ThreeVector(0.,HoleY_B,0.), FiberLV_B, "WLSFiber", ScintillatorLV, false, 1);

	auto FiberPV_A = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), FiberLV_A, "WLSFiber", iCladLV_A, false, 0);
	auto FiberPV_B = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), FiberLV_B, "WLSFiber", iCladLV_B, false, 1);

	//--------------------------------------------------
	// PhotonDet (Sensitive Detector)
	//--------------------------------------------------

	// Physical Construction
	auto SiPMS = new G4Tubs("SiPM", 0., WLSfiberR, SiPMThickness/2., 0.*deg, 360.*deg);

	auto SiPMLV_A0 = new G4LogicalVolume(SiPMS, FindMaterial("PMMA"), "SiPMA0LV");
	auto SiPMPV_A0 = new G4PVPlacement(0, G4ThreeVector(0., HoleY_A, (ScintZSize+SiPMThickness-CoatingThickness*2)/2), SiPMLV_A0, "SiPM", ExtrusionLV, false, 0);
	auto SiPMLV_A1 = new G4LogicalVolume(SiPMS, FindMaterial("PMMA"), "SiPMA1LV");
	auto SiPMPV_A1 = new G4PVPlacement(0, G4ThreeVector(0., HoleY_A, -(ScintZSize+SiPMThickness-CoatingThickness*2)/2), SiPMLV_A1, "SiPM", ExtrusionLV, false, 1);

	auto SiPMLV_B0 = new G4LogicalVolume(SiPMS, FindMaterial("PMMA"), "SiPMB0LV");
	auto SiPMPV_B0 = new G4PVPlacement(0, G4ThreeVector(0., HoleY_B, (ScintZSize+SiPMThickness-CoatingThickness*2)/2), SiPMLV_B0, "SiPM", ExtrusionLV, false, 2);
	auto SiPMLV_B1 = new G4LogicalVolume(SiPMS, FindMaterial("PMMA"), "SiPMB1LV");
	auto SiPMPV_B1 = new G4PVPlacement(0, G4ThreeVector(0., HoleY_B, -(ScintZSize+SiPMThickness-CoatingThickness*2)/2), SiPMLV_B1, "SiPM", ExtrusionLV, false, 3);

	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
	auto ScintillatorVA = new G4VisAttributes(G4Colour(0,0,0)); ScintillatorVA->SetVisibility(true); ScintillatorVA->SetForceWireframe(true);
	auto HoleVA= new G4VisAttributes(G4Colour(0,0,0,0.25)); HoleVA->SetVisibility(true); HoleVA->SetForceSolid(true);
	auto oCladVA= new G4VisAttributes(G4Colour(1,0.75,0.1,0.5)); oCladVA->SetVisibility(true); oCladVA->SetForceSolid(true);
	auto iCladVA= new G4VisAttributes(G4Colour(0.25,1.0,0.25,0.75)); iCladVA->SetVisibility(true); iCladVA->SetForceSolid(true);
	auto FiberVA= new G4VisAttributes(G4Colour(0.15,0.25,1.0,0.99)); FiberVA->SetVisibility(true); FiberVA->SetForceSolid(true);
	auto SiPMVA= new G4VisAttributes(G4Colour(1,0.1,0.25)); SiPMVA->SetVisibility(true); SiPMVA->SetForceSolid(true);

	// visual attributes for the shields 
	ExtrusionLV->SetVisAttributes(G4VisAttributes::GetInvisible()); 

	ScintillatorLV->SetVisAttributes(ScintillatorVA);
	HoleLV_A->SetVisAttributes(HoleVA);
	HoleLV_B->SetVisAttributes(HoleVA);
	oCladLV_A->SetVisAttributes(oCladVA);
	oCladLV_B->SetVisAttributes(oCladVA);
	iCladLV_A->SetVisAttributes(iCladVA);
	iCladLV_B->SetVisAttributes(iCladVA);
	FiberLV_A->SetVisAttributes(FiberVA);
	FiberLV_B->SetVisAttributes(FiberVA);
	SiPMLV_A0->SetVisAttributes(SiPMVA);
	SiPMLV_A1->SetVisAttributes(SiPMVA);
	SiPMLV_B0->SetVisAttributes(SiPMVA);
	SiPMLV_B1->SetVisAttributes(SiPMVA);

}


void DetectorConstruction::ConstructSDandField()
{
	G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

	std::string name[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};


	// declare scint as detector
	auto genDetector = new GenericSD("GenDet");
	G4SDManager::GetSDMpointer()->AddNewDetector(genDetector);
	SetSensitiveDetector("ScintillatorLV", genDetector);
	SetSensitiveDetector("FiberALV", genDetector);
	SetSensitiveDetector("FiberBLV", genDetector);
	auto scintDetector = new Scint_DetSD("ScintDet");
	G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);
	SetSensitiveDetector("SiPMA0LV", scintDetector);
	SetSensitiveDetector("SiPMA1LV", scintDetector);
	SetSensitiveDetector("SiPMB0LV", scintDetector);
	SetSensitiveDetector("SiPMB1LV", scintDetector);

}

G4Material* DetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}

void DetectorConstruction::DefineCommands()
{
	// Define /B5/generator command directory using generic messenger class
	fMessenger = new G4GenericMessenger(this, "/scint/", "Scintillator bar properties");

	// randomizePrimary command
	auto& lengthCmd = fMessenger->DeclarePropertyWithUnit("length", "cm", fLength);
	G4String guidance = "Scintillator bar length in cm.\n";
	lengthCmd.SetGuidance(guidance);
	lengthCmd.SetParameterName("length", true);
	lengthCmd.SetRange("length>=0.");
	lengthCmd.SetDefaultValue("200.");

	auto& widthCmd = fMessenger->DeclarePropertyWithUnit("width", "cm", fWidth);
	guidance = "Scintillator bar width in cm.\n";
	widthCmd.SetGuidance(guidance);
	widthCmd.SetParameterName("width", true);
	widthCmd.SetRange("width>=0.");
	widthCmd.SetDefaultValue("20.");

	auto& thicknessCmd = fMessenger->DeclarePropertyWithUnit("thickness", "cm", fThickness);
	guidance = "Scintillator bar thickness in cm.\n";
	thicknessCmd.SetGuidance(guidance);
	thicknessCmd.SetParameterName("thickness", true);
	thicknessCmd.SetRange("thickness>=0.");
	thicknessCmd.SetDefaultValue("2.");
	
	auto& yScaleCmd = fMessenger->DeclareProperty("yScale", fYScale);
	guidance = "Scaling factor for scintillator yield.\n";
	yScaleCmd.SetGuidance(guidance);
	yScaleCmd.SetParameterName("yScale", true);
	yScaleCmd.SetRange("yScale>=0.");
	yScaleCmd.SetDefaultValue("0.05");

}

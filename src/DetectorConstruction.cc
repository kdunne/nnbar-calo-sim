#include "DetectorConstruction.hh"
#include "detSD.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RotationMatrix.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4PSPopulation.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PhysicalConstants.hh"

#include <G4ProductionCuts.hh>
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include<string>
//#include "G4GDMLParser.hh"
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;



DetectorConstruction::DetectorConstruction()
	: G4VUserDetectorConstruction(), fCheckOverlaps(false)
{
}
//....
DetectorConstruction::~DetectorConstruction()
{}
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
	G4double a; G4double z; G4double density;

	G4Element* elH = nistManager->FindOrBuildElement("H");
	G4Element* elN = nistManager->FindOrBuildElement("N");
	G4Element* elO = nistManager->FindOrBuildElement("O");
	G4Element* elNa = nistManager->FindOrBuildElement("Na");
	G4Element* elAl = nistManager->FindOrBuildElement("Al");
	G4Element* elSi = nistManager->FindOrBuildElement("Si");
	G4Element* elP = nistManager->FindOrBuildElement("P");
	G4Element* elS = nistManager->FindOrBuildElement("S");
	G4Element* elK = nistManager->FindOrBuildElement("K");
	G4Element* elCa = nistManager->FindOrBuildElement("Ca");
	G4Element* elFe = nistManager->FindOrBuildElement("Fe");

	// Vacuum
	G4Material* vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);
	// -- optical properties of vacuum
	G4double vacuum_Energy[] = { 2.0 * eV,7.0 * eV, 7.14 * eV }; const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double); G4double vacuum_RIND[] = { 1.,1.,1. };
	G4MaterialPropertiesTable * vacuum_mt = new G4MaterialPropertiesTable();
	vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum); vacuum->SetMaterialPropertiesTable(vacuum_mt);

	G4Material* MagnadenseHC = new G4Material("MagnadenseHC",   density = 3.8*g/cm3,10);
	MagnadenseHC->AddElement(elH, 0.0053);
	MagnadenseHC->AddElement(elO, 0.332);
	MagnadenseHC->AddElement(elNa, 0.0046);
	MagnadenseHC->AddElement(elAl, 0.0064);
	MagnadenseHC->AddElement(elSi, 0.0469);
	MagnadenseHC->AddElement(elP, 0.0044);
	MagnadenseHC->AddElement(elS, 0.0003);
	MagnadenseHC->AddElement(elK, 0.0015);
	MagnadenseHC->AddElement(elCa, 0.0198);
	MagnadenseHC->AddElement(elFe, 0.579);

	// --------Air
	G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, 2);
	Air->AddElement(elN, 0.7); Air->AddElement(elO, 0.3);

}

//....
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	// Get materials
	auto defaultMaterial = G4Material::GetMaterial("Galactic");
	auto shieldMaterial = G4Material::GetMaterial("MagnadenseHC");

	if ( ! defaultMaterial || ! shieldMaterial) {
		G4ExceptionDescription msg; msg << "Cannot retrieve materials already defined.";
		G4Exception("DetectorConstruction::DefineVolumes()","MyCode0001", FatalException, msg);}

	// World
	auto worldSizeXY = 20.0 * m; auto worldSizeZ = 20.0 * m; //1 * calorThickness;
	auto worldS = new G4Box("WorldS",worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.);
	auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
	auto worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"WorldPV",0,false,0,fCheckOverlaps);

	G4double wall_half=7.*m;
	G4double det_half=5.*m;
	G4double inner_half=4.99*m;

	auto wall = new G4Box("wall",wall_half, wall_half, wall_half);
	auto det = new G4Box("det",det_half, det_half, det_half);
	auto inner = new G4Box("inner",inner_half, inner_half, inner_half);
		
	auto wallLV = new G4LogicalVolume(wall,shieldMaterial,"wallLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),wallLV,"wall",worldLV,false,0,fCheckOverlaps);  

	auto detLV = new G4LogicalVolume(det,defaultMaterial,"detLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),detLV,"det",wallLV,false,0,fCheckOverlaps);  

	auto innerLV = new G4LogicalVolume(inner,defaultMaterial,"innerLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),innerLV,"inner",detLV,false,0,fCheckOverlaps);  

	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
	auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
	auto green_color= new G4VisAttributes(G4Colour(0.517647,0.772549,0.556863)); green_color->SetVisibility(true);
	// original green: 
	auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
	auto red_color= new G4VisAttributes(G4Colour(0.956863,0.0901961,0.494118)); red_color->SetVisibility(true); 
	auto blue_color= new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765)); blue_color->SetVisibility(true);
	auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549)); black_color->SetVisibility(true);

	// visual attributes for the shields 
	return worldPV;
}

void DetectorConstruction::ConstructSDandField()
{ 

	G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
	G4String DetectorBoxName = "detLV" ;
	detSD* DetectorBox = new detSD(DetectorBoxName,0);
	G4SDManager::GetSDMpointer()->AddNewDetector(DetectorBox);
	SetSensitiveDetector("detLV", DetectorBox);
	
}

//....


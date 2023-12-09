#include "DetectorConstruction.hh"
#include "detSD.hh"
#include "samplingSD.hh"
#include "CVSD.hh"

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
	G4Element* elC = nistManager->FindOrBuildElement("C");
	G4Element* elTi = nistManager->FindOrBuildElement("Ti");
	G4Element* elAs = nistManager->FindOrBuildElement("As");
	G4Element* elPb = nistManager->FindOrBuildElement("Pb");
	G4Element* elO = nistManager->FindOrBuildElement("O");
	G4Element* elSi = nistManager->FindOrBuildElement("Si");
	G4Element* elN = nistManager->FindOrBuildElement("N");
	G4Element* elAr = nistManager->FindOrBuildElement("Ar");
	G4Element* elAl = nistManager->FindOrBuildElement("Al");
	G4Element* elMg = nistManager->FindOrBuildElement("Mg");
	G4Element* elMn = nistManager->FindOrBuildElement("Mn");
	G4Element* elNi = nistManager->FindOrBuildElement("Ni");
	G4Element* elCr = nistManager->FindOrBuildElement("Cr");
	G4Element* elF = nistManager->FindOrBuildElement("F");
	G4Element* elFe = nistManager->FindOrBuildElement("Fe");

	// Vacuum
	G4Material* vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);
	// -- optical properties of vacuum
	G4double vacuum_Energy[] = { 2.0 * eV,7.0 * eV, 7.14 * eV }; const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double); G4double vacuum_RIND[] = { 1.,1.,1. };
	G4MaterialPropertiesTable * vacuum_mt = new G4MaterialPropertiesTable();
	vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum); vacuum->SetMaterialPropertiesTable(vacuum_mt);

	// pipe
	G4Material* Al = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);

	G4Material* StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3,5);
	StainlessSteel->AddElement(elMn, 0.02);
	StainlessSteel->AddElement(elSi, 0.01);
	StainlessSteel->AddElement(elCr, 0.19);
	StainlessSteel->AddElement(elNi, 0.10);
	StainlessSteel->AddElement(elFe, 0.68);

	// Silicon
	G4Material* Silicon = new G4Material("Silicon", z=14., a= 28.0855*g/mole, density = 2.33*g/cm3);

	// --------Air
	G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, 2);
	Air->AddElement(elN, 0.7); Air->AddElement(elO, 0.3);

	// ---------FR4
	//----- Epoxy
	G4Material* Epoxy = new G4Material("Epoxy" , density=1.2*g/cm3, 2);
	Epoxy->AddElement(elH, 2); Epoxy->AddElement(elC, 2); 
	//----- SiO2 (Quarz)
	G4Material* SiO2 = new G4Material("SiO2",density= 2.200*g/cm3, 2);
	SiO2->AddElement(elSi, 1); SiO2->AddElement(elO , 2); 
	//FR4 (Glass + Epoxy)
	G4Material* FR4 = new G4Material("FR4" , density=1.86*g/cm3, 2);
	FR4->AddMaterial(Epoxy, 0.472); FR4->AddMaterial(SiO2, 0.528);

	//Carbon Target
	G4Material* Carbon_target = new G4Material("Carbon_target" , density=3.52*g/cm3, 1);
	Carbon_target->AddElement(elC, 1);

	// ----------TPC
	// CO2
	G4Material* CO2 = new G4Material("CO2", density= 1.98*g/cm3, 2);
	CO2->AddElement(elO, 2); CO2->AddElement(elC, 1);
	// Ar/CO2 80/20
	G4Material* Gas = new G4Material("Gas", density=1.823*mg/cm3, 2);
	Gas->AddElement(elAr, .8); Gas->AddMaterial(CO2, .2);

	// BC-408 taken from datasheet
	G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
	Scint->AddElement(elH, 0.524573); Scint->AddElement(elC, 1 - 0.524573);

	// Lead-glass (taken from PDG)
	G4Material* Abs = new G4Material("Abs", 6.22*g/cm3, 5);
	Abs->AddElement(elO, 0.156453); Abs->AddElement(elSi, 0.080866); Abs->AddElement(elTi, 0.008092); Abs->AddElement(elAs, .002651); Abs->AddElement(elPb, 0.751938);

	// MgF2 coating of the lead glass
	G4Material *AlMgF2 = new G4Material("AlMgF2", density = 2.9007*g/cm3, 3);
	AlMgF2->AddElement(elAl, 0.331); AlMgF2->AddElement(elF, 0.408); AlMgF2->AddElement(elMg, 0.261);
	G4Material* Lead = new G4Material("Lead", z=82., a= 207.2*g/mole, density = 11.29*g/cm3);


	// PMT window quartz actually
	G4Material* PMT_window = new G4Material("PMT_window_mat",density= 2.200*g/cm3, 2);
	PMT_window->AddElement(elSi, 1); PMT_window->AddElement(elO , 2); 
	// -- optical properties of PMT
	G4double pmt_window_Energy[] = { 2.0 * eV, 7.0 * eV}; const G4int pmt_window_num = sizeof(pmt_window_Energy) / sizeof(G4double); G4double pmt_window_RIND[] = { 1.53,1.53 };
	G4MaterialPropertiesTable * pmt_window_mt = new G4MaterialPropertiesTable();
	pmt_window_mt->AddProperty("RINDEX", pmt_window_Energy, pmt_window_RIND, pmt_window_num); PMT_window->SetMaterialPropertiesTable(pmt_window_mt);


	// ----------------- Generate and Add Material Properties Table ----------------
	G4double PhotonWavelength[] =
	{ 2325.4*nm, 1970.1*nm, 1529.6*nm, 1060.0*nm,
		1014.0*nm, 852.10*nm, 706.50*nm, 656.30*nm,
		643.80*nm, 632.80*nm, 589.30*nm, 587.60*nm,
		546.10*nm, 486.10*nm, 480.00*nm, 435.80*nm,
		404.70*nm, 365.00*nm};

	const G4int nEntries = sizeof(PhotonWavelength)/sizeof(G4double);
	G4double PhotonEnergy[nEntries];
	for (int i=0; i < nEntries; ++i) {PhotonEnergy[i] = (1240.*nm/PhotonWavelength[i])*eV;};

	//// Lead Glass Schott SF5
	G4double refractiveIndex[] =
	{ 1.63289, 1.63785, 1.64359, 1.65104,
		1.65206, 1.65664, 1.66327, 1.66661,
		1.66756, 1.66846, 1.67252, 1.67270,
		1.67764, 1.68750, 1.68876, 1.69986,
		1.71069, 1.73056};

	G4MaterialPropertiesTable* absMPT = new G4MaterialPropertiesTable();
	absMPT->AddProperty("RINDEX", PhotonEnergy, refractiveIndex, nEntries,false,true);
	Abs->SetMaterialPropertiesTable(absMPT);
	//G4cout << "Absorber Properties -------" << G4endl;
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	absMPT->DumpTable();

	//Scintillator Optical Properties

	const G4int nEntries2 = 12;

	G4double ScintPhotonEnergy[nEntries2] =

	{ 2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV, 2.76*eV,
		2.82*eV, 2.92*eV, 2.95*eV, 3.02*eV, 3.1*eV,
		3.26*eV, 3.44*eV};

	G4double rindex_scint[nEntries2] =
	{1.58, 1.58, 1.58, 1.58, 1.58,
		1.58, 1.58, 1.58, 1.58, 1.58,
		1.58, 1.58};

	G4double atten_scint[nEntries2] =
	{210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
		210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
		210*cm, 210*cm};

	G4double scintilFast[nEntries2] = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};
	G4double scintilSlow[nEntries2] = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};

	G4MaterialPropertiesTable *scintMPT = new G4MaterialPropertiesTable();
	scintMPT->AddProperty("RINDEX", ScintPhotonEnergy, rindex_scint, nEntries2,false,true);
	scintMPT->AddProperty("ABSLENGTH", ScintPhotonEnergy, atten_scint, nEntries2,false,true);
	scintMPT->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhotonEnergy, scintilFast, nEntries2,false,true);
	scintMPT->AddProperty("SCINTILLATIONCOMPONENT2", ScintPhotonEnergy, scintilSlow, nEntries2,false,true);
	// 64% of Antracene: 17400
	scintMPT->AddConstProperty("SCINTILLATIONYIELD", 11136./ MeV); //original 11136000.
	scintMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
	scintMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 0.9*ns); // org: 0.9
	scintMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 2.1*ns); // org: 2.1
	scintMPT->AddConstProperty("SCINTILLATIONYIELD1", 0.5);
	scintMPT->AddConstProperty("SCINTILLATIONYIELD2", 0.5);
	Scint->SetMaterialPropertiesTable(scintMPT);
	scintMPT->DumpTable();
}

//....
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	// Get materials
	auto defaultMaterial = G4Material::GetMaterial("Galactic");
	auto scintMaterial = G4Material::GetMaterial("Scint");
	auto shieldMaterial = G4Material::GetMaterial("StainlessSteel");
	auto attenuatorMaterial = G4Material::GetMaterial("Galactic");

	if ( ! defaultMaterial || ! scintMaterial) {
		G4ExceptionDescription msg; msg << "Cannot retrieve materials already defined.";
		G4Exception("DetectorConstruction::DefineVolumes()","MyCode0001", FatalException, msg);}

	// World
	auto worldSizeXY = 10.0 * m; auto worldSizeZ = 10.0 * m; //1 * calorThickness;
	auto worldS = new G4Box("WorldS",worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.);
	auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
	auto worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"WorldPV",0,false,0,fCheckOverlaps);

// HIBEAM
	G4double bar_l = 4.2*m;
	G4double bar_w = 3.2*m;
	G4double bar_t = 2.0*cm;
	G4double shortbar_l = 1.2*m;
// NNBAR
//	G4double bar_l = 6.4*m;
//	G4double bar_t = 3.0*cm;
//	G4double shortbar_l = 2.1*m;

	G4double detector_half=bar_w/2.+2*bar_t+1*mm;
	G4double detector_length=bar_l/2.+2*bar_t+1*mm;
	auto detector = new G4Box("detector",detector_half, detector_half, detector_length);
	auto detectorLV = new G4LogicalVolume(detector,defaultMaterial,"detectorLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),detectorLV,"detector",worldLV,false,0,fCheckOverlaps);  

	// Build rudimentary annihilation detector
		
// HIBEAM
	G4double passive_half=1.45*m;
	G4double passive_length=2.1*m;
	G4double sampling_half=1.26*m;
	G4double sampling_length=2.00*m;
	G4double lg_half=1.21*m;
	G4double lg_length=1.95*m;
	G4double scint_half=0.905*m;
	G4double scint_length=1.65*m;
	G4double tpc_half=0.555*m;
	G4double tpc_length=1.27*m;
	G4double pipe_outer_radius=0.395*m;
	G4double pipe_inner_radius=0.375*m;
	//G4double pipe_length=6.0*m;
	
	auto passive_box = new G4Box("passive",passive_half, passive_half, passive_length);
	auto sampling_box = new G4Box("sampling",sampling_half, sampling_half, sampling_length);
	auto leadglass_box = new G4Box("leadglass",lg_half, lg_half, lg_length);
	auto scintillator_box = new G4Box("scintillator",scint_half, scint_half, scint_length);
	auto tpc_box = new G4Box("scintillator", tpc_half, tpc_half, tpc_length);
	auto gap_box = new G4Box("gap", pipe_outer_radius+1*mm, pipe_outer_radius+1*mm, passive_length);
		
	auto passive0 = new G4SubtractionSolid("passive0",passive_box, sampling_box);
	auto passive = new G4SubtractionSolid("passive",passive0,gap_box);
	auto passiveLV = new G4LogicalVolume(passive,shieldMaterial,"passiveLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),passiveLV,"passive",detectorLV,false,0,fCheckOverlaps);  

	auto sampling = new G4SubtractionSolid("sampling",sampling_box,leadglass_box);
	auto samplingLV = new G4LogicalVolume(sampling,defaultMaterial,"samplingLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),samplingLV,"sampling",detectorLV,false,0,fCheckOverlaps);  

	auto leadglass0 = new G4SubtractionSolid("leadglass0",leadglass_box,scintillator_box);
	auto leadglass = new G4SubtractionSolid("leadglass",leadglass0,gap_box);
	auto leadglassLV = new G4LogicalVolume(leadglass,G4Material::GetMaterial("Abs"),"leadglassLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),leadglassLV,"leadglass",detectorLV,false,0,fCheckOverlaps);  

	auto scintillator0 = new G4SubtractionSolid("scintillator",scintillator_box,tpc_box);
	auto scintillator = new G4SubtractionSolid("scintillator",scintillator0,gap_box);
	auto scintillatorLV = new G4LogicalVolume(scintillator,scintMaterial,"scintillatorLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),scintillatorLV,"scintillator",detectorLV,false,0,fCheckOverlaps);  

	auto tpc = new G4SubtractionSolid("tpc",tpc_box,gap_box);
	auto tpcLV = new G4LogicalVolume(tpc,G4Material::GetMaterial("Gas"),"tpcLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tpcLV,"tpc",detectorLV,false,0,fCheckOverlaps);  
		
	auto beampipe = new G4Cons("beampipe", pipe_inner_radius, pipe_outer_radius, pipe_inner_radius, pipe_outer_radius,sampling_length,0.,360.*deg);
	auto beampipeLV = new G4LogicalVolume(beampipe,G4Material::GetMaterial("Aluminum"),"beampipeLV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),beampipeLV,"beampipeLV",detectorLV,false,0,fCheckOverlaps);  

	auto beampipe2 = new G4Cons("beampipe2", pipe_inner_radius, pipe_outer_radius, pipe_inner_radius, pipe_outer_radius,(detector_length-sampling_length)/2.,0.,360.*deg);
	auto beampipe2LV = new G4LogicalVolume(beampipe2,G4Material::GetMaterial("Aluminum"),"beampipe2LV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,sampling_length+(detector_length-sampling_length)/2.),beampipe2LV,"beampipe2LV",detectorLV,false,0,fCheckOverlaps);  

	auto beampipe3 = new G4Cons("beampipe3", pipe_inner_radius, pipe_outer_radius, pipe_inner_radius, pipe_outer_radius,(detector_length-sampling_length)/2.,0.,360.*deg);
	auto beampipe3LV = new G4LogicalVolume(beampipe2,G4Material::GetMaterial("Aluminum"),"beampipe3LV");
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-sampling_length-(detector_length-sampling_length)/2.),beampipe3LV,"beampipe3LV",detectorLV,false,0,fCheckOverlaps);  
	// ***********************************************
	// ======= Scintillator CV building  =========
	// ***********************************************
	G4double pipe_gap = pipe_outer_radius + 0.02*m;
	G4double attenuator_gap = 1*cm;
	// Scint CV up/down
	G4double veto_ud_x = bar_w; 
	G4double veto_ud_y = bar_t-0.05*cm; 
	G4double veto_ud_z = bar_l;

	auto CV_ud_1_S = new G4Box("CV_ud_S",veto_ud_x/2., veto_ud_y/2., veto_ud_z/2.);
	auto CV_ud_1_LV = new G4LogicalVolume(CV_ud_1_S,defaultMaterial,"CV_ud_1_LV");

	auto CV_ud_2_S = new G4Box("CV_ud_S",veto_ud_x/2., veto_ud_y/2., veto_ud_z/2.);
	auto CV_ud_2_LV = new G4LogicalVolume(CV_ud_2_S,defaultMaterial,"CV_ud_2_LV");

	auto Att_ud_S = new G4Box("Att_ud_S",veto_ud_x/2., attenuator_gap/2., veto_ud_z/2.);
	auto Att_ud_LV = new G4LogicalVolume(Att_ud_S,attenuatorMaterial,"Att_ud_LV");

	new G4PVPlacement(0,G4ThreeVector(0.,bar_w/2.+bar_t/2.,0.),CV_ud_1_LV,"up_layer1",detectorLV,false,0,fCheckOverlaps);  
	new G4PVPlacement(0,G4ThreeVector(0.,bar_w/2.+3.*bar_t/2.+attenuator_gap,0.),CV_ud_2_LV,"up_layer2",detectorLV,false,1,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(0.,-(bar_w/2.+bar_t/2.),0.),CV_ud_1_LV,"down_layer1",detectorLV,false,2,fCheckOverlaps);  
	new G4PVPlacement(0,G4ThreeVector(0.,-(bar_w/2.+3.*bar_t/2.+attenuator_gap),0.),CV_ud_2_LV,"down_layer2",detectorLV,false,3,fCheckOverlaps);
	
	new G4PVPlacement(0,G4ThreeVector(0.,bar_w/2.+bar_t+attenuator_gap/2.,0.),Att_ud_LV,"up_attenuator",detectorLV,false,101,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(0.,-(bar_w/2.+bar_t+attenuator_gap/2.),0.),Att_ud_LV,"down_attenuator",detectorLV,false,102,fCheckOverlaps);

	// ===> Filling Scint CV ud with scint bars
	int num_ud_xbars = 16;
	int num_ud_zbars = 21;

	G4double xbar_ud_w = veto_ud_x/num_ud_xbars*mm;
	G4double zbar_ud_w = veto_ud_z/num_ud_zbars*mm;

	auto CV_xbar_ud_S = new G4Box("CV_ud_S",xbar_ud_w/2., bar_t/2., bar_l/2.);
	auto CV_xbar_ud_LV = new G4LogicalVolume(CV_xbar_ud_S,scintMaterial,"CV_xbar_ud_LV");

	for (int i = 0;i<num_ud_xbars;i++){
		new G4PVPlacement(0,G4ThreeVector(-veto_ud_x/2.+((2.0*i+1)*xbar_ud_w/2.),0.,0),CV_xbar_ud_LV,"xbar"+std::to_string(i),CV_ud_1_LV,false,i,fCheckOverlaps);
	}

	auto CV_zbar_ud_S = new G4Box("CV_ud_S",bar_w/2., bar_t/2., zbar_ud_w/2.);
	auto CV_zbar_ud_LV = new G4LogicalVolume(CV_zbar_ud_S,scintMaterial,"CV_zbar_ud_LV");

	for (int i = 0;i<num_ud_zbars;i++){
		new G4PVPlacement(0,G4ThreeVector(0.,0.,-veto_ud_z/2.+((2.0*i+1)*zbar_ud_w/2.)),CV_zbar_ud_LV,"xbar"+std::to_string(i),CV_ud_2_LV,false,i,fCheckOverlaps);
	}

	// Scint CV lr
	G4double veto_lr_x = bar_t-0.05*cm; 
	G4double veto_lr_y = bar_w; 
	G4double veto_lr_z = bar_l;

	auto CV_lr_1_S = new G4Box("CV_lr_S",veto_lr_x/2., veto_lr_y/2., veto_lr_z/2.);
	auto CV_lr_1_LV = new G4LogicalVolume(CV_lr_1_S,defaultMaterial,"CV_lr_1_LV");

	auto CV_lr_2_S = new G4Box("CV_lr_S",veto_lr_x/2., veto_lr_y/2., veto_lr_z/2.);
	auto CV_lr_2_LV = new G4LogicalVolume(CV_lr_2_S,defaultMaterial,"CV_lr_2_LV");

	auto Att_lr_S = new G4Box("Att_lr_S",attenuator_gap/2., veto_lr_y/2., veto_lr_z/2.);
	auto Att_lr_LV = new G4LogicalVolume(Att_lr_S,attenuatorMaterial,"Att_lr_LV");

	new G4PVPlacement(0,G4ThreeVector(bar_w/2.+bar_t/2.,0.,0.),CV_lr_1_LV,"left_layer1",detectorLV,false,4,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(bar_w/2.+3.*bar_t/2.+attenuator_gap,0.,0.),CV_lr_2_LV,"left_layer2",detectorLV,false,5,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(-(bar_w/2.+bar_t/2.),0.,0.),CV_lr_1_LV,"right_layer1",detectorLV,false,6,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(-(bar_w/2.+3.*bar_t/2.+attenuator_gap),0.,0.),CV_lr_2_LV,"right_layer2",detectorLV,false,7,fCheckOverlaps);

	new G4PVPlacement(0,G4ThreeVector(bar_w/2.+bar_t+attenuator_gap,0.,0.),Att_lr_LV,"left_attenuator",detectorLV,false,103,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(-(bar_w/2.+bar_t+attenuator_gap),0.,0.),Att_lr_LV,"right_attenuator",detectorLV,false,104,fCheckOverlaps);
	
	// ===> Filling Scint CV top with scint bars
	int num_lr_ybars = 16;
	int num_lr_zbars = 21;

	G4double ybar_lr_w = veto_lr_y/num_lr_ybars*mm;
	G4double zbar_lr_w = veto_lr_z/num_lr_zbars*mm;

	auto CV_ybar_lr_S = new G4Box("CV_lr_S",bar_t/2.,ybar_lr_w/2., bar_l/2.);
	auto CV_ybar_lr_LV = new G4LogicalVolume(CV_ybar_lr_S,scintMaterial,"CV_ybar_lr_LV");

	for (int i = 0;i<num_lr_ybars;i++){
		new G4PVPlacement(0,G4ThreeVector(0.,-veto_lr_y/2.+((2.0*i+1)*ybar_lr_w/2.),0.),CV_ybar_lr_LV,"xbar"+std::to_string(i),CV_lr_1_LV,false,i,fCheckOverlaps);
	}

	auto CV_zbar_lr_S = new G4Box("CV_top_S",bar_t/2., bar_w/2., zbar_lr_w/2.);
	auto CV_zbar_lr_LV = new G4LogicalVolume(CV_zbar_lr_S,scintMaterial,"CV_zbar_lr_LV");

	for (int i = 0;i<num_lr_zbars;i++){
		new G4PVPlacement(0,G4ThreeVector(0.,0.,-veto_lr_z/2.+((2.0*i+1)*zbar_lr_w/2.)),CV_zbar_lr_LV,"xbar"+std::to_string(i),CV_lr_2_LV,false,i,fCheckOverlaps);
	}

	// Scint CV front/back
	G4double veto_fb_x = bar_w; 
	G4double veto_fb_y = bar_w; 
	G4double veto_fb_z = bar_t-0.05*cm;

	auto CV_fb_1_S = new G4Box("CV_fb_S",veto_fb_x/2., veto_fb_y/2., veto_fb_z/2.);
	auto CV_fb_1_LV = new G4LogicalVolume(CV_fb_1_S,defaultMaterial,"CV_fb_1_LV");

	auto CV_fb_2_S = new G4Box("CV_fb_S",veto_fb_x/2., veto_fb_y/2., veto_fb_z/2.);
	auto CV_fb_2_LV = new G4LogicalVolume(CV_fb_2_S,defaultMaterial,"CV_fb_2_LV");

	auto Att_fb_0 = new G4Box("Att_fb_0",veto_fb_x/2., veto_fb_y/2., veto_fb_z/2.);
	auto Att_fb_S = new G4SubtractionSolid("Att_fb_S",Att_fb_0,gap_box);
	auto Att_fb_LV = new G4LogicalVolume(Att_fb_S,attenuatorMaterial,"Att_fb_LV");

	new G4PVPlacement(0,G4ThreeVector(0.,0.,bar_l/2.+bar_t/2.),CV_fb_1_LV,"back_layer1",detectorLV,false,8,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(0.,0.,bar_l/2.+3.*bar_t/2.+attenuator_gap),CV_fb_2_LV,"back_layer2",detectorLV,false,9,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-(bar_l/2.+bar_t/2.)),CV_fb_1_LV,"front_layer1",detectorLV,false,10,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-(bar_l/2.+3.*bar_t/2.+attenuator_gap)),CV_fb_2_LV,"front_layer2",detectorLV,false,11,fCheckOverlaps);

	new G4PVPlacement(0,G4ThreeVector(0.,0.,bar_l/2.+bar_t+attenuator_gap/2.),Att_fb_LV,"back_attenuator",detectorLV,false,105,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-(bar_l/2.+bar_t+attenuator_gap/2.)),Att_fb_LV,"front_attenuator",detectorLV,false,106,fCheckOverlaps);
	// ===> Filling Scint CV top with scint bars
	int num_fb_xbars = 16;
	int num_fb_ybars = 16;

	G4double xbar_fb_w = veto_fb_x/num_fb_xbars*mm;
	G4double ybar_fb_w = veto_fb_y/num_fb_ybars*mm;

	auto CV_xbar_fb_S = new G4Box("CV_fb_S",xbar_fb_w/2.,bar_w/2., bar_t/2.);
	auto CV_xbar_fb_LV = new G4LogicalVolume(CV_xbar_fb_S,scintMaterial,"CV_xbar_fb_LV");
	
	auto CV_xshortbar_fb_S = new G4Box("CV_fb_S",xbar_fb_w/2.,shortbar_l/2., bar_t/2.);
	auto CV_xshortbar_fb_LV = new G4LogicalVolume(CV_xshortbar_fb_S,scintMaterial,"CV_xshortbar_fb_LV");

	for (int i = 0;i<num_fb_xbars;i++){
		G4double pos = -veto_fb_x/2.+((2.0*i+1)*xbar_fb_w/2.);
		std::cout << "xbar pos: " << pos << std::endl;
		G4double edge[2] = {pos-xbar_fb_w/2.,pos+xbar_fb_w/2.};
		if(edge[1]<-pipe_gap || edge[0]>pipe_gap){
			new G4PVPlacement(0,G4ThreeVector(pos,0.,0.),CV_xbar_fb_LV,"xbar"+std::to_string(i),CV_fb_1_LV,false,i,fCheckOverlaps);
		}
		else{
			new G4PVPlacement(0,G4ThreeVector(pos,-0.5*(bar_w-shortbar_l),0.),CV_xshortbar_fb_LV,"xbar"+std::to_string(i)+"a",CV_fb_1_LV,false,i,fCheckOverlaps);
			new G4PVPlacement(0,G4ThreeVector(pos,0.5*(bar_w-shortbar_l),0.),CV_xshortbar_fb_LV,"xbar"+std::to_string(i)+"b",CV_fb_1_LV,false,i,fCheckOverlaps);
		}
	}

	auto CV_ybar_fb_S = new G4Box("CV_top_S",bar_w/2., ybar_fb_w/2., bar_t/2.);
	auto CV_ybar_fb_LV = new G4LogicalVolume(CV_ybar_fb_S,scintMaterial,"CV_ybar_fb_LV");

	auto CV_yshortbar_fb_S = new G4Box("CV_top_S",shortbar_l/2., ybar_fb_w/2., bar_t/2.);
	auto CV_yshortbar_fb_LV = new G4LogicalVolume(CV_yshortbar_fb_S,scintMaterial,"CV_yshortbar_fb_LV");
	
	for (int i = 0;i<num_fb_ybars;i++){
		G4double pos = -veto_fb_y/2.+((2.0*i+1)*ybar_fb_w/2.);
		std::cout << "ybar pos: " << pos << std::endl;
		G4double edge[2] = {pos-ybar_fb_w/2.,pos+ybar_fb_w/2.};
		if(edge[1]<-pipe_gap || edge[0]>pipe_gap){
			new G4PVPlacement(0,G4ThreeVector(0.,pos,0.),CV_ybar_fb_LV,"ybar"+std::to_string(i),CV_fb_2_LV,false,i,fCheckOverlaps);
		}
		else{
			new G4PVPlacement(0,G4ThreeVector(-0.5*(bar_w-shortbar_l),pos,0.),CV_yshortbar_fb_LV,"ybar"+std::to_string(i)+"a",CV_fb_2_LV,false,i,fCheckOverlaps);
			new G4PVPlacement(0,G4ThreeVector(0.5*(bar_w-shortbar_l),pos,0.),CV_yshortbar_fb_LV,"ybar"+std::to_string(i)+"b",CV_fb_2_LV,false,i,fCheckOverlaps);
		}
	}

	// cuts for shield
	G4Region* shield_region = new G4Region("CV_region");

	shield_region->AddRootLogicalVolume(CV_xbar_ud_LV);
	shield_region->AddRootLogicalVolume(CV_zbar_ud_LV);
	shield_region->AddRootLogicalVolume(CV_ybar_lr_LV);
	shield_region->AddRootLogicalVolume(CV_zbar_lr_LV);
	shield_region->AddRootLogicalVolume(CV_xbar_fb_LV);
	shield_region->AddRootLogicalVolume(CV_ybar_fb_LV);

	G4Region* shield_Region = G4RegionStore::GetInstance()->GetRegion("CV_region");
	G4ProductionCuts* shieldcut = new G4ProductionCuts();
	shieldcut->SetProductionCut(8.0*cm,"gamma"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
	shieldcut->SetProductionCut(5.0*mm,"e-");
	shieldcut->SetProductionCut(5.0*mm,"e+");
	shieldcut->SetProductionCut(15.0*mm,"proton");
	shield_Region->SetProductionCuts(shieldcut);

	//new G4LogicalSkinSurface("name",CV_bar_top_1_LV,op_glass_world);
	//new G4LogicalSkinSurface("name",worldLV,op_glass_world);


	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
	auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true); grey_color->SetForceSolid(true);
	auto green_color= new G4VisAttributes(G4Colour(0.517647,0.772549,0.556863)); green_color->SetVisibility(true);
	// original green: 
	auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
	auto red_color= new G4VisAttributes(G4Colour(0.956863,0.0901961,0.494118)); red_color->SetVisibility(true); red_color->SetForceWireframe(true); 
	auto blue_color= new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765,0.5)); blue_color->SetVisibility(true); blue_color->SetForceSolid(true);
	auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549)); black_color->SetVisibility(true); black_color->SetForceWireframe(true); black_color->SetLineWidth(2);

	// visual attributes for the shields 
	leadglassLV->SetVisAttributes(blue_color);
	scintillatorLV->SetVisAttributes(black_color);
	tpcLV->SetVisAttributes(black_color);
	beampipeLV->SetVisAttributes(grey_color);
	beampipe2LV->SetVisAttributes(grey_color);
	beampipe3LV->SetVisAttributes(grey_color);
//	leadglassLV->SetVisAttributes(black_color);
//	scintillatorLV->SetVisAttributes(black_color);
//	tpcLV->SetVisAttributes(black_color);
//	beampipeLV->SetVisAttributes(black_color);
	CV_ud_1_LV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
	CV_ud_2_LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CV_lr_1_LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CV_lr_2_LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CV_fb_1_LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CV_fb_2_LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	Att_ud_LV->SetVisAttributes(black_color); 
	Att_lr_LV->SetVisAttributes(black_color); 
	Att_fb_LV->SetVisAttributes(black_color); 
	
	CV_xbar_ud_LV->SetVisAttributes(red_color);
	CV_zbar_ud_LV->SetVisAttributes(red_color);
	CV_ybar_lr_LV->SetVisAttributes(red_color);
	CV_zbar_lr_LV->SetVisAttributes(red_color);
	CV_xbar_fb_LV->SetVisAttributes(red_color);
	CV_ybar_fb_LV->SetVisAttributes(red_color);
	CV_xshortbar_fb_LV->SetVisAttributes(red_color);
	CV_yshortbar_fb_LV->SetVisAttributes(red_color);

	return worldPV;
}

//....

void DetectorConstruction::ConstructSDandField()
{ 

	G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
	G4String DetectorBoxName = "detectorLV" ;
	detSD* DetectorBox = new detSD(DetectorBoxName,1);
	G4SDManager::GetSDMpointer()->AddNewDetector(DetectorBox);
	SetSensitiveDetector("detectorLV", DetectorBox);
	
	G4String samplingBoxName = "samplingLV" ;
	samplingSD* samplingBox = new samplingSD(samplingBoxName);
	G4SDManager::GetSDMpointer()->AddNewDetector(samplingBox);
	SetSensitiveDetector("samplingLV", samplingBox);

	G4String CVDetectorName = "CVLV" ;
	CVSD* shieldDetector = new CVSD(CVDetectorName);
	G4SDManager::GetSDMpointer()->AddNewDetector(shieldDetector);
	SetSensitiveDetector("CV_xbar_ud_LV", shieldDetector);
	SetSensitiveDetector("CV_zbar_ud_LV", shieldDetector);
	SetSensitiveDetector("CV_ybar_lr_LV", shieldDetector);
	SetSensitiveDetector("CV_zbar_lr_LV", shieldDetector);
	SetSensitiveDetector("CV_xbar_fb_LV", shieldDetector);
	SetSensitiveDetector("CV_ybar_fb_LV", shieldDetector);

}

//....


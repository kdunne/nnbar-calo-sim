#include "WLSMaterials.hh"
#include "G4SystemOfUnits.hh"

WLSMaterials* WLSMaterials::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSMaterials::WLSMaterials()
{
  fNistMan = G4NistManager::Instance();
  fNistMan->SetVerbose(2);
  CreateMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSMaterials::~WLSMaterials()
{
  delete    fPMMA;
  delete    fPethylene;
  delete    fFPethylene;
  delete    fPolystyrene;
  delete    fSilicone;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSMaterials* WLSMaterials::GetInstance()
{
  if (fInstance == 0)
    {
      fInstance = new WLSMaterials();
    }
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* WLSMaterials::GetMaterial(const G4String material)
{
  G4Material* mat =  fNistMan->FindOrBuildMaterial(material);

  if (!mat) mat = G4Material::GetMaterial(material);
  if (!mat) {
     std::ostringstream o;
     o << "Material " << material << " not found!";
     G4Exception("WLSMaterials::GetMaterial","",
                 FatalException,o.str().c_str());
  } else {
    std::cout << "Material Found: " << material << std::endl;
  }

  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSMaterials::CreateMaterials()
{
  G4double density;
  G4int ncomponents;
  G4double fractionmass;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  // Materials Definitions
  // =====================

  //--------------------------------------------------
  // Vacuum
  //--------------------------------------------------

  fNistMan->FindOrBuildMaterial("G4_Galactic");

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  fAir = fNistMan->FindOrBuildMaterial("G4_AIR");
 
  //--------------------------------------------------
  // Si
  //--------------------------------------------------
 
  fSi = fNistMan->FindOrBuildMaterial("G4_Si");

  //--------------------------------------------------
  // WLSfiber PMMA
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  density = 1.190*g/cm3;

  fPMMA = fNistMan->
          ConstructNewMaterial("PMMA", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Cladding (polyethylene)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.200*g/cm3;

  fPethylene = fNistMan->
          ConstructNewMaterial("Pethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Double Cladding (fluorinated polyethylene)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.400*g/cm3;

  fFPethylene = fNistMan->
          ConstructNewMaterial("FPethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Polystyrene
  //--------------------------------------------------
 
  elements.push_back("C");     natoms.push_back(8);
  elements.push_back("H");     natoms.push_back(8);

  density = 1.050*g/cm3;

  fPolystyrene = fNistMan->
          ConstructNewMaterial("Polystyrene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Silicone (Template for Optical Grease)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(6);
  
  density = 1.060*g/cm3;

  fSilicone = fNistMan->
          ConstructNewMaterial("Silicone", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Aluminium
  //--------------------------------------------------

  fNistMan->FindOrBuildMaterial("G4_Al");

  //--------------------------------------------------
  // TiO2
  //--------------------------------------------------

  elements.push_back("Ti");     natoms.push_back(1);
  elements.push_back("O");      natoms.push_back(2);

  density     = 4.26*g/cm3;

  G4Material* TiO2 = fNistMan->
          ConstructNewMaterial("TiO2", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Scintillator Coating - 15% TiO2 and 85% polystyrene by weight.
  //--------------------------------------------------

  density = 1.52*g/cm3;

  fCoating =
          new G4Material("Coating", density, ncomponents=2);

  fCoating->AddMaterial(TiO2,         fractionmass = 15*perCent);
  fCoating->AddMaterial(fPolystyrene, fractionmass = 85*perCent);

  //
  // ------------ Generate & Add Material Properties Table ------------
  //

  G4double photonEnergy[] =
  {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
   2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
   2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
   2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
   2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
   2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
   2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
   3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
   3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
   3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  G4double refractiveIndex[] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
  
  assert(sizeof(refractiveIndex) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);

  fAir->SetMaterialPropertiesTable(mpt);


  //--------------------------------------------------
  // Si
  //--------------------------------------------------

  G4double refractiveIndexSi[] =
  { 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88,
    3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88,
    3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88,
    3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88,
    3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88, 3.88};

  G4double absLengthSi[] =
  { 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm,
    0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm,
    0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm,
    0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm,
    0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm, 0.01*mm};

  assert(sizeof(refractiveIndex) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* mptSi = new G4MaterialPropertiesTable();
  mptSi->AddProperty("RINDEX", photonEnergy, refractiveIndexSi, nEntries);
  mptSi->AddProperty("ABSLENGTH", photonEnergy, absLengthSi, nEntries);

  fSi->SetMaterialPropertiesTable(mptSi);

  //--------------------------------------------------
  //  PMMA for WLSfibers
  //--------------------------------------------------

  std::vector<G4double> wls_Energy = {2.00*eV, 2.16*eV, 2.25*eV, 2.36*eV, 2.43*eV, 2.48*eV, 2.53*eV, 2.61*eV, 2.76*eV, 2.92*eV, 3.1*eV, 3.5*eV};
  std::vector<G4double> refractiveIndexWLSfiber =  { 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59};
  std::vector<G4double> absFiber = {  3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m, 3.0*m };
  std::vector<G4double> WLSabsFiber = { 10.0*m, 10.0*m, 10.0*m, 10.0*m, 10.0*m, 2.0*mm, 0.1*mm, 0.1*mm, 0.1*mm, 1.0*mm, 5.0*mm, 10.0*m };
  std::vector<G4double> emissionFib = { 0.0, 0.3, 0.6, 0.7, 1., 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  
  // Add entries into properties table
  G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
  mptWLSfiber->AddProperty("RINDEX",wls_Energy,refractiveIndexWLSfiber);
  mptWLSfiber->AddProperty("ABSLENGTH",wls_Energy,absFiber);
  mptWLSfiber->AddProperty("WLSABSLENGTH",wls_Energy,WLSabsFiber);
  mptWLSfiber->AddProperty("WLSCOMPONENT",wls_Energy,emissionFib);
  mptWLSfiber->AddConstProperty("WLSTIMECONSTANT", 2.5*ns);

  fPMMA->SetMaterialPropertiesTable(mptWLSfiber);

  //--------------------------------------------------
  //  Polyethylene
  //--------------------------------------------------

  std::vector<G4double> wls_Energy_short = {2.00*eV, 3.5*eV};
  std::vector<G4double> refractiveIndexClad1 = { 1.49, 1.49,};
  std::vector<G4double> absClad = {20.0*m,20.0*m};

  // Add entries into properties table
  G4MaterialPropertiesTable* mptClad1 = new G4MaterialPropertiesTable();
  mptClad1->AddProperty("RINDEX",wls_Energy_short,refractiveIndexClad1);
  mptClad1->AddProperty("ABSLENGTH",wls_Energy_short,absClad);

  fPethylene->SetMaterialPropertiesTable(mptClad1);

  //--------------------------------------------------
  // Fluorinated Polyethylene
  //--------------------------------------------------

  std::vector<G4double> refractiveIndexClad2 = { 1.42, 1.42};

  // Add entries into properties table
  G4MaterialPropertiesTable* mptClad2 = new G4MaterialPropertiesTable();
  mptClad2->AddProperty("RINDEX",wls_Energy_short,refractiveIndexClad2);
  mptClad2->AddProperty("ABSLENGTH",wls_Energy_short,absClad);

  fFPethylene->SetMaterialPropertiesTable(mptClad2);

  //--------------------------------------------------
  //  Polystyrene
  //--------------------------------------------------

  std::vector<G4double> ScintPhotonEnergy =  { 2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV, 2.76*eV, 2.82*eV, 2.92*eV, 2.95*eV, 3.02*eV, 3.1*eV, 3.26*eV, 3.44*eV};
  std::vector<G4double> rindex_scint = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58};
  std::vector<G4double> atten_scint = {210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm};
  std::vector<G4double> scintilFast = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};
  std::vector<G4double> scintilSlow = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};

  G4MaterialPropertiesTable *mptPolystyrene = new G4MaterialPropertiesTable();
  mptPolystyrene->AddProperty("RINDEX", ScintPhotonEnergy, rindex_scint,false,true);
  mptPolystyrene->AddProperty("ABSLENGTH", ScintPhotonEnergy, atten_scint,false,true);
  mptPolystyrene->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhotonEnergy, scintilFast,false,true);
  mptPolystyrene->AddProperty("SCINTILLATIONCOMPONENT2", ScintPhotonEnergy, scintilSlow,false,true);
  // 64% of Antracene: 17400
  double yield = 17400*0.64;
  double yieldScale = 0.0125;
  mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD", yield*yieldScale/MeV); //original 11136.
  mptPolystyrene->AddConstProperty("RESOLUTIONSCALE", sqrt(yieldScale));
  mptPolystyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 0.9*ns); // org: 0.9
  mptPolystyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 2.1*ns); // org: 2.1
  mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD1", 1.);
  fPolystyrene->SetMaterialPropertiesTable(mptPolystyrene);
  
  fPolystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4cout << "Scintillator Properties -------" << G4endl;
  mptPolystyrene->DumpTable();
}

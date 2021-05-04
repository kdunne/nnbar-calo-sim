#include "PrimaryGeneratorAction.hh"
#include <iomanip>
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "Analysis.hh"
//.....

using namespace std;

int particle_name_file_index;
extern std::ofstream Particle_outFile;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern G4double event_number;
extern G4int run_number;

PrimaryGeneratorAction::PrimaryGeneratorAction()
	:fParticleGun(nullptr)
{
	const char* inputfile = "setup.file";
	fMessenger = new G4GenericMessenger(this, "/particle_generator/", "Name the particle for the file name");
	G4GenericMessenger::Command& filenameCMD = fMessenger->DeclareProperty("Particle_index", particle_name_file_index, "Index of the particle in the order");
	filenameCMD.SetParameterName("Particle Index", true);
	filenameCMD.SetDefaultValue("99");

	fParticleGun = new G4ParticleGun();

	std::ifstream inputFile;
	inputFile.open(inputfile, std::ios::in);

	char buffer[1000];

	if (inputFile.fail()) {
		if (*inputfile != 0)  //....only complain if a filename was given
			std::cerr << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << std::endl;
		InputState = -1;
	}
	else {
		std::string setupString("");
		while (!inputFile.getline(buffer, 1000).eof()) {
			setupString.append(buffer);
			setupString.append(" ");
		}

		CRYSetup* setup = new CRYSetup(setupString, "/home/billy/nnbar/cry_v1.7/data");

		gen = new CRYGenerator(setup);

		// set random number generator
		RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
		setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
		InputState = 0;

	}
	// create a vector to store the CRY particle properties
	vect = new std::vector<CRYParticle*>;
	particleTable = G4ParticleTable::GetParticleTable();
	gunMessenger = new PrimaryGeneratorMessenger(this);
	std::cerr << "Input state: " << InputState << std::endl;
	std::cerr << "particle table: " << particleTable << std::endl;
}

//....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{}

//....

void PrimaryGeneratorAction::InputCRY()
{
	InputState = 1;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::UpdateCRY(std::string * MessInput)
{
	CRYSetup* setup = new CRYSetup(*MessInput, "/home/billy/nnbar/cry_v1.7/data");

	gen = new CRYGenerator(setup);

	// set random number generator
	RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
	setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
	InputState = 0;

}

void PrimaryGeneratorAction::CRYFromFile(G4String newValue)
{
	// Read the cry input file
	std::ifstream inputFile;
	inputFile.open(newValue, std::ios::in);
	char buffer[1000];

	if (inputFile.fail()) {
		std::cerr << "Failed to open input file " << newValue << std::endl;
		std::cerr << "Make sure to define the cry library on the command line" << std::endl;
		InputState = -1;
	}
	else {

		std::string setupString("");
		while (!inputFile.getline(buffer, 1000).eof()) {
			setupString.append(buffer);
			setupString.append(" ");
		}

		CRYSetup* setup = new CRYSetup(setupString, "/home/billy/nnbar/cry_v1.7/data");

		gen = new CRYGenerator(setup);

		// set random number generator
		RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
		setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
		InputState = 0;
	}

	std::cerr << "Input state after CRYFromFile: " << InputState << std::endl;

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{

	if (InputState != 0) {
		G4String* str = new G4String("CRY library was not successfully initialized");
		std::cerr << "Error in Generate Primaries" << std::endl;
		G4Exception("PrimaryGeneratorAction", "1", RunMustBeAborted, *str);
	}

	G4double x; G4double y; G4double z;
	G4double t; G4double px; G4double py; G4double pz;
	G4double KE;

	G4String particleName;
	vect->clear();
	gen->genEvent(vect);

	for (unsigned j = 0; j < vect->size(); j++) {
		

		std::vector<G4double> particle_gun_record_row;

		particleName = CRYUtils::partName((*vect)[j]->id());

		x = (*vect)[j]->x()* m;
		y = 5.1*m;
		z = (*vect)[j]->y()* m; // (*vect)[j]->z() * m
		KE = (*vect)[j]->ke() * MeV;
		px = (*vect)[j]->u();
		py = (*vect)[j]->v();
		pz = (*vect)[j]->w();
		t = (*vect)[j]->t(); 

		G4String all_type[7] = { "neutron","proton","gamma","electron","muon","pion","kaon" };
		G4double name_ID = 99;

		for (int i = 0; i <= 6; i++) { if (particleName == all_type[i]) { name_ID = i; break; } }

		fParticleGun->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
		fParticleGun->SetParticleEnergy(KE);
		fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz)); //(*vect)[j]->w())
		fParticleGun->SetParticleTime((*vect)[j]->t());
		
		std::cout << event_number << "  " << particleName << " ID: " << (*vect)[j]->PDGid() << " charge= " << (*vect)[j]->charge() << " "
			<< setprecision(4)
			<< " energy (MeV)=" << (*vect)[j]->ke() * MeV << " "
			<< " pos (m)"
			<< G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), z / m)
			<< " " << "direction cosines "
			<< G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w())
			<< " " << "Particle Time: " << (*vect)[j]->t() 
		<< std::endl;


		Particle_outFile <<  event_number << ",";
		Particle_outFile <<  (*vect)[j]->PDGid() << ","; //PID
		Particle_outFile << particleTable->FindParticle((*vect)[j]->PDGid())-> GetPDGMass() << ",";
		Particle_outFile <<  0 << ","; // charge 
		Particle_outFile <<  KE << ","; // 100+n*10 
		Particle_outFile <<  x << ",";
		Particle_outFile <<  y << ",";
		Particle_outFile <<  z << ",";
		Particle_outFile <<  t << ",";
		Particle_outFile <<  px << ",";
		Particle_outFile <<  py << ",";
		Particle_outFile <<  pz << G4endl;

		fParticleGun->GeneratePrimaryVertex(anEvent);
		
		delete (*vect)[j];
	
	}
}

//....


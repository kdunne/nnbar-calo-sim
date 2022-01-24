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
{fParticleGun = new G4ParticleGun();}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
	G4double x; G4double y; G4double z;
	G4double t; G4double px; G4double py; G4double pz;
	G4double KE;

	G4String particleName;

	std::vector<G4double> particle_gun_record_row;

	int n = run_number; // every energy, it has 50 entries


	double radius = 0.0;
	double x_rand = radius*G4UniformRand();

	x = x_rand*m;
	y = sqrt(pow(radius,2.0)-pow(x_rand,2.0))*m; //y = 0.0*m;
	z = 0.0* m; // (*vect)[j]->z() * m
	KE = 50.0+n*25.0  *MeV; //250.0 * MeV; +n*50 
	px = G4UniformRand();//G4UniformRand();
	py = G4UniformRand();//G4UniformRand();
	pz = G4UniformRand();//G4UniformRand();
	t = 0.; 

	G4String all_type[7] = { "neutron","proton","gamma","electron","muon","pion","kaon" };
	G4double name_ID = 99;

	for (int i = 0; i <= 6; i++) { if (particleName == all_type[i]) { name_ID = i; break; } }

	auto particle_name = "proton";
	auto particle_ID = 2212;

	fParticleGun->SetParticleDefinition(particleTable->FindParticle(particle_name));
	fParticleGun->SetParticleEnergy(KE);	
	fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz)); //(*vect)[j]->w())
	fParticleGun->SetParticleTime(0.0);
	Particle_outFile <<  event_number << ",";
	Particle_outFile <<  particle_ID << ","; //PID
	Particle_outFile <<  particleTable -> FindParticle(particle_name) -> GetPDGMass() << ","; //Mass
	Particle_outFile <<  particleTable -> FindParticle(particle_name) -> GetPDGCharge() << ","; // charge 
	Particle_outFile <<  50+n*25 << ","; // +n*50 
	Particle_outFile <<  x << ",";
	Particle_outFile <<  y << ",";
	Particle_outFile <<  z << ",";
	Particle_outFile <<  t << ",";
	Particle_outFile <<  px << ",";
	Particle_outFile <<  py << ",";
	Particle_outFile <<  pz << G4endl;

	fParticleGun->GeneratePrimaryVertex(anEvent);

	std::cout<< "Event Number -- " << event_number << "position: (" << x << "," << y << "," << z 
	<< ") Direction: (" << px << "," << py << "," << pz << ")"
	<< std::endl;

}

//....


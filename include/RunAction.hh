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

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <string>
#include <vector>

using namespace std;

class G4Run;
extern std::vector<int> event_ID;
extern std::vector<int> particle_ID; // particle type (represented by a number)
extern std::vector<double> particle_KE; //initial KE of the particle
extern std::vector<double> particle_momentum_x; // Initial momentum (direction info included)
extern std::vector<double> particle_momentum_y;
extern std::vector<double> particle_momentum_z;
extern std::vector<double> particle_time; // Time stamp in second
extern std::vector<double> particle_x; // initial position of the particle
extern std::vector<double> particle_y;
extern std::vector<double> particle_z;
extern std::vector<std::vector<int>> scint_photon_collection_pri;
extern std::vector<std::vector<int>> scint_photon_collection_all;
extern std::vector<int> lead_glass_photon_all;
extern std::vector<int> lead_glass_photon_pri;


class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

  private:
	std::vector<string> particle_name{ "Neutron","Proton","Gamma","Electron","Muon","Pion","Kaon" };



    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
};

//....

#endif


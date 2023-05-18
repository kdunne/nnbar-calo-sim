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

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "G4THitsMap.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

class HistoManager;

class EventAction : public G4UserEventAction
{
	public:
		EventAction(HistoManager *histo);
		virtual ~EventAction();

		virtual void  BeginOfEventAction(const G4Event* event);
		virtual void    EndOfEventAction(const G4Event* event);

	private:
		HistoManager* fHistoManager;

		const G4int nPlanes=12;
		const G4int nBarsPerPlane=64;
		const G4int nBars=nPlanes*nBarsPerPlane;
		const G4double bar_length = 640*cm;
		const G4double bar_thickness = 3*cm;
		const G4double bar_width = 10*cm;

		// methods
		G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
				const G4Event* event) const;
	//	G4double CalcEnergy(G4double edep, G4double w, G4double r);
	//	G4double CalcTime(G4double w, G4double v);
		void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
				G4double gapEdep, G4double scintTrackLength) const; //, G4double gapTrackLength) const;

		// data members                   
		G4int  detHitsCollectionID;
		G4int  samplingHitsCollectionID;
		G4int  CVHitsCollectionID;

	};

//....

#endif



#ifndef BarHit_h
#define BarHit_h 1

#include "globals.hh"
#include "CLHEP/Units/PhysicalConstants.h"

class BarHit
{
public:
    BarHit();
    BarHit(G4int p, G4int b, G4double d, G4double l, G4double w, G4double t);
    BarHit(const BarHit&);
    virtual ~BarHit();
 
    const BarHit& operator=(const BarHit&);
    G4bool operator==(const BarHit&) const;

private:
	
	const G4double c_scint=CLHEP::c_light/1.59; //from Kuraray data sheet
	const G4double AttLen=3.0*m; //from Kuraray data sheet
		
	const G4double ppMeV=17400*0.56;
	const G4double parr[3] = {0.139, 110.68*mm, 0.0687};
	const G4double parw[2] = {0.2395, 1860.*mm};

	G4int plane;
	G4int bar;
	
	G4double distance;	// distance of surface from center of target
	G4double length;	
	G4double width;
	G4double thickness;
    
	G4double x0;
	G4double y0;

    G4double time;
    G4double x;
    G4double y;
    G4double z;
    G4double eDep;

	G4double e1, e2, e3, e4;
	G4double t1, t2, t3, t4;
	G4double post, pose;
    
public:
    void SetGeom(G4int p, G4int b, G4double d, G4double l, G4double w, G4double t);
	void AddHit(G4double xx, G4double yy, G4double zz, G4double ee);
	void AnalyzeHits();
	
	G4double CalcEnergy(G4double ed, G4double w, G4double r);
	G4double CalcTime(G4double w, G4double v);

    inline G4int GetPlane(){return plane;}
    inline G4int GetBar(){return bar;}
    
	inline G4double GetTime(){return time;}
    inline G4double GetX(){return x;}
    inline G4double GetY(){return y;}
    inline G4double GetZ(){return z;}
    inline G4double GetEDep(){return eDep;}
    	
	inline G4double GetE1(){return e1;}
    inline G4double GetE2(){return e2;}
    inline G4double GetE3(){return e3;}
    inline G4double GetE4(){return e4;}

	inline G4double GetT1(){return t1;}
    inline G4double GetT2(){return t2;}
    inline G4double GetT3(){return t3;}
    inline G4double GetT4(){return t4;}

    inline G4double GetPosT(){return post;}
    inline G4double GetPosE(){return pose;}
    
	inline void SetTime(G4double t){time = t;}
   // inline void SetX(G4double x0){x = x0;}
   // inline void SetY(G4double y0){y = y0;}
   // inline void SetZ(G4double z0){z = z0;}
   // inline void SetEDep(G4double e){eDep = e;}

};

#endif

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "CLHEP/Random/JamesRandom.h"

#include "BarHit.hh"

BarHit::BarHit()
{
	plane = 0;
	bar = 0;
	distance = 50*cm;
	length = 100*cm;
	width = 10*cm;
	thickness = 2*cm;
    time = 0.;
	x = 0.;
	y = 0.;
	z = 0.;
    eDep = 0.;	

	x0 = distance + thickness*(plane+0.5);
	y0 = distance - width*(bar+0.5);

}

BarHit::BarHit(G4int p, G4int b, G4double d, G4double l, G4double w, G4double t)
{
	plane = p;
	bar = b;
	distance = d;
	length = l;
	width = w;
	thickness = t;
    time = 0.;
	x = 0.;
	y = 0.;
	z = 0.;
    eDep = 0.;
	
	x0 = distance + thickness*(plane+0.5);
	y0 = distance - width*(bar+0.5);
}
BarHit::~BarHit()
{}

BarHit::BarHit(const BarHit& right)
{
    time = right.time;
    x = right.x;
    y = right.y;
    z = right.z;
    eDep = right.eDep;
}

const BarHit& BarHit::operator=(const BarHit& right)
{
    time = right.time;
	x = right.x;
	y = right.y;
    z = right.z;
	eDep = right.eDep;

    return *this;
}

G4bool BarHit::operator==(const BarHit& right) const
{
	return((time==right.time)&&(x==right.x)&&(y==right.y)&&(z==right.z)&&(eDep==right.eDep));
}

void BarHit::SetGeom(G4int p, G4int b, G4double d, G4double l, G4double w, G4double t)
{
	plane = p;
	bar = b;
	distance = d;
	length = l;
	width = w;
	thickness = t;
	
	x0 = distance + thickness*(plane%2+0.5);
	y0 = distance - width*(bar+0.5);
}

void BarHit::AddHit(G4double xx, G4double yy, G4double zz, G4double ee)
{
	eDep += ee;

	xx = (xx>0) ? xx-x0 : xx+x0;
	yy = yy+y0;
	x += (xx*ee); 
	y += (yy*ee); 
	z += (zz*ee);
}


G4double BarHit::CalcEnergy(G4double ed, G4double w, G4double r)
{
	
	CLHEP::HepRandomEngine* engine = new CLHEP::HepJamesRandom();
	CLHEP::RandGaussQ rGauss(engine);
  	
	r = (r<10*cm) ? r : 10*cm;
	G4double np = floor(ed*ppMeV);
	np *= (parr[0]*exp(-r/parr[1])+parr[2]);
	np *= (parw[0]*exp(-w/parw[1]));
	np = rGauss.shoot(np,sqrt(np));
	np = (np<0) ? -np : np;

	return np/ppMeV;  
}  

G4double BarHit::CalcTime(G4double w, G4double v)
{
	
	CLHEP::HepRandomEngine* engine = new CLHEP::HepJamesRandom();
	CLHEP::RandGaussQ rGauss(engine);
  	
	G4double t=sqrt(w*w+v*v)/c_scint;
	t = rGauss.shoot(t,0.5*ns);
	t = (t<0) ? -t : t;
	return t;  
}  

void BarHit::AnalyzeHits()
{
	x /= eDep;
	y /= eDep;
	z /= eDep;

	G4double z1 = z+length/2;
	G4double z2 = length/2-z;
	G4double y1 = y+width/4.;
	G4double y2 = y-width/4.;
	G4double r1 = sqrt(x*x+y1*y1);
	G4double r2 = sqrt(x*x+y2*y2);

	t1 = CalcTime(z1,y1);
	t2 = CalcTime(z2,y1);
	t3 = CalcTime(z1,y2);
	t4 = CalcTime(z2,y2);

	e1 = CalcEnergy(eDep, z1, r1);
	e2 = CalcEnergy(eDep, z2, r1);
	e3 = CalcEnergy(eDep, z1, r2);
	e4 = CalcEnergy(eDep, z2, r2);

	post = 0.5*(t1+t3-t2-t4)*c_scint;
	pose = 0.5*log((e2+e4)/(e1+e3))*AttLen;
}


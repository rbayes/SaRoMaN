//
// **************************************************************************
// * Mind test toroidal field                                               *
// *                                                                        *
// * This is a class based on the ExN04Field class. The purpose is to test  *
// * the impacts of using a toroidal magnetic field within the Mind         *
// * simulation. For simulation of a more realistic detector a field map    *
// * should be measured and or simulated and sampled as needed.             *
// *                                                  Ryan Bayes, 12.04.11  *
// **************************************************************************

#include "babyMindField.hh"
#include <CLHEP/Units/GlobalSystemOfUnits.h>

babyMindField::babyMindField(double Bmag, double height, 
			     double width, int npanels)
{
  Bx = Bmag * tesla;
  hmax = height*cm;
  wmax = width*cm;
  seg = npanels;
}

babyMindField::~babyMindField()
{;}

void babyMindField::GetFieldValue(const double Point[3], double *Bfield) const 
{
  // use the components of the Bfield vector to set radial, azimuthal,
  // and longitudinal components within a cylindrical geometry
  double x = Point[0];
  double y = Point[1];
  double z = Point[2];
  double r = sqrt(x*x + y*y);

  // Only define the magnetic field within the detector and ignore it
  // outside of the detector. 
  
  if( x > hmax/2. && y > hmax/2. ) {
    Bfield[0] = 0.;
    Bfield[1] = 0.;
    Bfield[2] = 0.;
  } else {
    if ( seg == 2 ){
      // the field is divided into two equal sections
      if ( y > 0 ){
	Bfield[0] = Bx;
	Bfield[1] = 0.;
	Bfield[2] = 0.;
      } else { 
	Bfield[0] = -Bx;
	Bfield[1] = 0.;
	Bfield[2] = 0.;
      }
    } else if ( seg > 2 ){
      // There are two outside plates half the size of the inner plates
      double w = hmax/2./double(seg - 1);
      int iseg = int(y / w);
      if ( iseg % 2 == 0 ){
	Bfield[0] = Bx;
	Bfield[1] = 0.;
	Bfield[2] = 0.;
      } else { 
	Bfield[0] = -Bx;
	Bfield[1] = 0.;
	Bfield[2] = 0.;
      }
    }
  }
}

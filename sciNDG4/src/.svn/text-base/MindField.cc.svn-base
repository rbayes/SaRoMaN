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

#include "MindField.hh"
#include <CLHEP/Units/GlobalSystemOfUnits.h>

MindField::MindField(double Bmag, double height, double width, bool doublew)
{
  Btheta = Bmag * tesla;
  hmax = height*cm;
  rmin =   50*cm;
  wmax = width*cm;
  isLBNO = doublew;
}

MindField::~MindField()
{;}

void MindField::GetFieldValue(const double Point[3], double *Bfield) const 
{
  // use the components of the Bfield vector to set radial, azimuthal,
  // and longitudinal components within a cylindrical geometry
  double x = Point[0];
  double y = Point[1];
  double z = Point[2];
  double r = sqrt(x*x + y*y);

  // Only define the magnetic field within the detector and ignore it
  // outside of the detector. 
  
  double xoct = hmax * ( 1. + 1./(1. + sqrt(2.)))/2.;
  if(isLBNO) xoct += wmax/4.; 
  if(!isLBNO){
    // azimuthal field Bx = -Bmag*sin(theta), By = Bmag*cos(theta) The
    // radial dependence is cooked to look a little like the field map
    // along the x axis 
    Bfield[0] = -Btheta * ( 1.36 + 0.0406*m/r + 0.80*exp(-r*0.16/m)) * y/r;
    Bfield[1] =  Btheta * ( 1.36 + 0.0406*m/r + 0.80*exp(-r*0.16/m)) * x/r;
    Bfield[2] =  0.;
  }
  else{
    double rL = sqrt(pow(x-wmax/4,2) + y*y);
    double rR = sqrt(pow(x+wmax/4,2) + y*y); 
    Bfield[0] = -Btheta * ( 1.36 + 0.0406*m/rL + 0.80*exp(-rL*0.16/m)) * y/rL
       + Btheta * ( 1.36 + 0.0406*m/rR + 0.80*exp(-rR*0.16/m)) * y/rR;
    Bfield[1] =  Btheta * ( 1.36 + 0.0406*m/rL + 0.80*exp(-rL*0.16/m)) * (x-wmax/4.)/rL
      - Btheta * ( 1.36 + 0.0406*m/rR + 0.80*exp(-rR*0.16/m)) * (x+wmax/4.)/rR;
    Bfield[2] =  0.;
  }

}

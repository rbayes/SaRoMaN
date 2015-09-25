
#ifndef MindField_H
#define MindField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class MindField : public G4MagneticField
{
public:
  MindField(double Bmag, double height, double length, bool doublew=false);
  ~MindField();

  void GetFieldValue( const double Point[3], double *Bfield) const;

private:
  double Btheta;
  double hmax;
  double rmin;
  double wmax;
  bool isLBNO;
};

#endif

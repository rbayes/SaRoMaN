
#ifndef babyMindField_H
#define babyMindField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class babyMindField : public G4MagneticField
{
public:
  babyMindField(double Bmag, double height, double length, int npanels);
  ~babyMindField();

  void GetFieldValue( const double Point[3], double *Bfield) const;

private:
  double Bx;
  double hmax;
  double wmax;
  int seg;
};

#endif

#ifndef MindFieldMapR_H
#define MindFieldMapR_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include <string>
#include <vector>
#include <map>
#include <math.h>

class MindFieldMapR : public G4MagneticField{
protected:
  double* _xmax;
  double* _xmin;
  double* _ymax;
  double* _ymin;
public:
  MindFieldMapR(G4String Bmap, double fieldScaling, double passive_thickness,
		double number_active, double active_thickness);
  ~MindFieldMapR();
  
  void GetFieldValue(const double Point[3], double *Bfield) const; 
  //void SortMapValues(const double Point[3], double *Bfield); 
  
private:
  void Tokenize(const std::string& str,
		std::vector<std::string>& tokens,
		const std::string& delimiters=" ");

  double _gridLimit;
  bool gridtested;
  double _fieldScale;
  double _passive_thickness;
  double _number_active;
  double _active_thickness;
  double dx;
  double dy;
  double Nx;
  double Ny;
  double _Xmin;
  double _Xmax;
  double _Ymin;
  double _Ymax;
  std::map<double, std::map<double, std::vector<double> > > dtmap;
  std::vector<std::vector<std::vector<double> > > Bvec;
  std::vector<double> xvec;
  std::vector<double> yvec;
};

#endif

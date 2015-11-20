//
// **************************************************************************
// * Mind Field Map Reader
// *                                                                        *
// * This is a class based on the ExN04Field class. The underlying          *
// * is that the magnetic field is an ordered list of six numbers           *
// * [x,y,z,Bx,By,Bz]. There should be some form of linear, or parabolic    *
// * interpolation between points (depends on how fast the field changes).  *
// *                                                  Ryan Bayes, 28.04.11  *
// **************************************************************************

#include "MindFieldMapR.hh"
#include "G4String.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//#include <CLHEP/Units/GlobalSystemOfUnits.h>
#include "G4SystemOfUnits.hh"

using namespace std;

MindFieldMapR::MindFieldMapR(G4String Bmap, double fieldScaling, 
			     double passive_thickness, double number_active,
			     double active_thickness)
  :_fieldScale(fieldScaling), _gridLimit(8*cm), gridtested(false),
   _passive_thickness(passive_thickness), _number_active(number_active),
   _active_thickness(active_thickness)
{
  
  ifstream bfile;
  bfile.open(Bmap.c_str(), ifstream::in);
  
  char temp[256];
  _Xmin = 999999999.9*m;
  _Xmax =-999999999.9*m;
  _Ymin = 999999999.9*m;
  _Ymax =-999999999.9*m;

  while( bfile.getline( temp, 256 ) ){
    // Fill a vector with data points defined by temp
    std::string s(temp);
    std::vector<std::string> entries;
    Tokenize(s, entries);
    // if(entries.size() != 6) G4err
    double x = atof(entries.at(0).c_str())*m;
    double y = atof(entries[1].c_str())*m;
    double z = atof(entries[2].c_str())*m;
    double Bx = atof(entries[3].c_str())*tesla;
    double By = atof(entries[4].c_str())*tesla;
    double Bz = atof(entries[5].c_str())*tesla;
    std::vector<double> B;
    B.push_back(Bx);
    B.push_back(By);
    B.push_back(Bz);
    dtmap[x][y] = B;
    if(x < _Xmin){ _Xmin = x; }
    if(x > _Xmax){ _Xmax = x; }
    if(y < _Ymin){ _Ymin = y; }
    if(y > _Ymax){ _Ymax = y; }
  }
  Nx = dtmap.size();
  Ny = dtmap[_Xmax].size();
  dx = (_Xmax - _Xmin)/dtmap.size();
  dy = (_Ymax - _Ymin)/dtmap[_Xmax].size();
  
  Bvec.reserve(dtmap.size());
  std::map<double, std::map<double, std::vector<double> > >::iterator ix;
  std::map<double, std::vector<double> >::iterator iy;
  for( ix=dtmap.begin(); ix != dtmap.end(); ix++ ){
    xvec.push_back((*ix).first);
    std::vector< std::vector<double> > ytemp;
    // std::cout<<(*ix).first<<"\t"<<(*ix).second.size()<<std::endl;
    for( iy=(*ix).second.begin(); iy != (*ix).second.end(); iy++){
      if(ix == dtmap.begin()) yvec.push_back((*iy).first);
      ytemp.push_back((*iy).second);
    }
    Bvec.push_back(ytemp);
  }
  // Print statement to tell us if the field is actually read.
  std::cout<<"Field Map Consists of "<< dtmap.size()
	   <<" data points, for "<<_Xmin<<" < x < "<<_Xmax
	   <<" and "<<dtmap[_Xmin].size()<<" data points, for "
	   <<_Ymin<<" < y < "<<_Ymax<< std::endl;
  // std::cout<<tesla<<std::endl;
  bfile.close();
}

MindFieldMapR::~MindFieldMapR()
{;}

void MindFieldMapR::GetFieldValue(const double Point[3], double *Bfield) const 
{
  //  double piece_length = _passive_thickness + _number_active*_active_thickness;
  //  double z = fmod((Point[2] + piece_length/2.), piece_length);
  //  bool inFe = z < _passive_thickness ? true: false;

  if((Point[0] < _Xmin + dx || Point[0] > _Xmax - dx 
      || Point[1] < _Ymin + dy || Point[1] > _Ymax - dy)){
    Bfield[0] = 0;
    Bfield[1] = 0;
    Bfield[2] = 0;
  }
  else{
    // Assume a regular grid of points (could this be a problem????)
    int j1 = int(floor((Point[0] - _Xmin)/dx + 0.5));
    int k1 = int(floor((Point[1] - _Ymin)/dy + 0.5));
    if(j1 < 0) j1 = 0;
    if(k1 < 0) k1 = 0;
    if(j1 >= Bvec.size())     j1 = Bvec.size() - 1;
    if(k1 >= Bvec[j1].size()) k1 = Bvec[j1].size() - 1;
    double x1 = _Xmin + double(j1)*dx;
    double x2 = _Xmin + double(j1-1)*dx;
    double y1 = _Ymin + double(k1)*dy;
    double y2 = _Ymin + double(k1-1)*dy;
    int j2 = j1-1;
    int k2 = k1-1;
    if ((Point[0] > x1 && j1 <= Nx-2) || j1 == 0){
      x2 = _Xmin + double(j1 + 1)*dx;
      j2 = j1+1;
    }
    if ((Point[1] > y1 && k1 <= Ny-2) || k1 == 0){
      y2 = _Ymin + double(k1 + 1)*dy;
      k2 = k1+1;
    }
    
    double Bx0 = Bvec[j1][k1][0];
    double By0 = Bvec[j1][k1][1];
    double Bx1 = Bvec[j2][k1][0];
    double By1 = Bvec[j2][k1][1];
    /*if( (Bx1 < 0.5*Bx0 || By1 < 0.5*By0) && j2 < j1 ){ 
      // A correction to avoid the edges of the plate
      j2 = j1+1;
      x2 = _Xmin + double(j2)*dx;
      Bx1 = Bvec[j2][k1][0];
      By1 = Bvec[j2][k1][1];
    }
    else if( (Bx1 < 0.5*Bx0 || By1 < 0.5*By0) && j2 > j1 ){
      j2 = j1-1;
      x2 = _Xmin + double(j2)*dx;
      Bx1 = Bvec[j2][k1][0];
      By1 = Bvec[j2][k1][1];
      }*/
    double Bx2 = Bvec[j1][k2][0];
    double By2 = Bvec[j1][k2][1];
    /*if( (Bx2 < 0.5*Bx0 || By2 < 0.5*By0) && k2 < k1 ){ 
      // A correction to avoid the edges of the plate
      k2 = k1+1;
      y2 = _Ymin + double(k2)*dx;
      Bx1 = Bvec[j1][k2][0];
      By1 = Bvec[j1][k2][1];
    }
    else if( (Bx2 < 0.5*Bx0 || By2 < 0.5*By0) && k2 > k1 ){
      k2 = k1-1;
      y2 = _Ymin + double(k2)*dx;
      Bx1 = Bvec[j1][k2][0];
      By1 = Bvec[j1][k2][1];
      }*/
    
    double dBxdx = (Bx1 - Bx0)/(x2 - x1);
    double dBydx = (By1 - By0)/(x2 - x1);
    double dBxdy = (Bx2 - Bx0)/(y2 - y1);
    double dBydy = (By2 - By0)/(y2 - y1);
    // if (fabs(x1 - Point[0]) < fabs(x2 - Point[0]) ) ix = j1;
    // if (fabs(y1 - Point[1]) < fabs(y2 - Point[1]) ) iy = k1;
    // Interpolation in x, y, and z removed as it was not sufficient near slots.
    Bfield[0] = _fieldScale*(Bx0 + dBxdx*(Point[0] - x1) + dBxdy*(Point[1] - y1));
    Bfield[1] = _fieldScale*(By0 + dBydx*(Point[0] - x1) + dBydy*(Point[1] - y1));
    Bfield[2] = _fieldScale*(Bvec[j1][k1][2]);
    
    //std::cout<<Point[0]<<"\t"<<Point[1]<<"\t"<<Point[2]
    // <<Bfield[0]<<"\t"<<Bfield[1]<<"\t"<<Bfield[2]<<std::endl;
  }
}

void MindFieldMapR::Tokenize(const std::string& str,
			     std::vector<std::string>& tokens,
			     const std::string& delimiters)
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos){
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}





/*
void MindFieldMapR::SortMapValues(const double Point[3], double *Bfield){
  // Reinitiallize reference fieldpoints
  MindFieldPoint minX(Point[0],Point[1],Point[2],0.0,0.0,0.0,_gridLimit);
  MindFieldPoint maxX(Point[0],Point[1],Point[2],0.0,0.0,0.0,_gridLimit);
  MindFieldPoint minY(Point[0],Point[1],Point[2],0.0,0.0,0.0,_gridLimit);
  MindFieldPoint maxY(Point[0],Point[1],Point[2],0.0,0.0,0.0,_gridLimit);
  //The following is a failed attempt to make the algorithm more efficient
  //  following the observation that a particle does not necessarily move 
  //  by  more than the 5cm file grid spacing. Will not work because geant4
  //  requires the function to be declaired constant and a constant function
  //  cannot alter members of the class
  double cdgXmin =  9999999.9;
  double cdgXmax =  9999999.9;
  double cdgYmin =  9999999.9;
  double cdgYmax =  9999999.9;
  if(gridtested){
    // Check distance to reference points to determine if the loop
    // is necessary.
    cdgXmin = sqrt(pow(Point[0] - _xmin[0],2) +
		   pow(Point[1] - _xmin[1],2));
    cdgXmax = sqrt(pow(Point[0] - _xmax[0],2) +
		   pow(Point[1] - _xmax[1],2));
    cdgYmin = sqrt(pow(Point[0] - _ymin[0],2) +
		   pow(Point[1] - _ymin[1],2));
    cdgYmax = sqrt(pow(Point[0] - _ymax[0],2) +
		   pow(Point[1] - _ymax[1],2));
  }
  if(cdgXmin < _gridLimit && cdgXmax < _gridLimit && 
     cdgYmin < _gridLimit && cdgYmax < _gridLimit){
    minX.Update(_xmin);
    maxX.Update(_xmax);
    minY.Update(_xmin);
    maxY.Update(_ymax);
  }
  else {
    
    // Collect all the data points within 8 cm of the test point and
    // interpolate the magnetic field for the closest point Grid is
    // probably small enough and the field is smooth enough that linear
    // interpolation between the nearest neighbours is likely sufficient
    
    // Could not use iterators to streamline/optimize function because of 
    // required const function definition
    bool x0=false, x1=false, y0=false, y1=false;
    for(int k = 0; k < points.size(); k++){
      bool used = false;
      //while(it != points.end()){
      double cdg = sqrt(pow(Point[0] - points[k].xx(),2) +
			pow(Point[1] - points[k].yy(),2));
      if(cdg < _gridLimit){
	if(minX.xx() > points[k].xx() && !x0 && !used) {
	  minX = points[k]; x0 = true; used = true; continue; }
	if(maxX.xx() < points[k].xx() && !x1 && !used) {
	  maxX = points[k]; x1 = true; used = true; continue; }
	if(minY.yy() > points[k].yy() && !y0 && !used) {
	  minY = points[k]; y0 = true; used = true; continue; }
	if(maxY.yy() < points[k].yy() && !y1 && !used) {
	  maxY = points[k]; y1 = true; used = true; continue; }
	// If I get this far I should break from the loop
	if(x0 && x1 && y0 && y1) break;
      }
    }
  }
  
  
  //  std::cout<<minX.xx()<<"\t"<<maxX.xx()<<"\t"
  //  <<minX.yy()<<"\t"<<maxX.yy()<<std::endl;
  //  
  //  std::cout<<minY.xx()<<"\t"<<maxY.xx()<<"\t"
  //  <<minY.yy()<<"\t"<<maxY.yy()<<std::endl;
  
  //Now to conduct the interpolation with respect to x-y plane
  double dbxx = (maxX.Bx() - minX.Bx())/(maxX.xx() - minX.xx());
  double dbxy = (maxY.Bx() - minY.Bx())/(maxY.yy() - minY.yy());
  double dbyx = (maxX.By() - minX.By())/(maxX.xx() - minX.xx());
  double dbyy = (maxY.By() - minY.By())/(maxY.yy() - minY.yy());
  // Planer interpolation of magnetic field for Bx
  Bfield[0] = (dbxx * (Point[0] - minX.xx()) + 
	       dbxy * (Point[1] - minX.yy()) + 
	       minX.Bx());
  // Planer interpolation of magnetic field for By
  Bfield[1] = (dbyx * (Point[0] - minX.xx()) + 
	       dbyy * (Point[1] - minX.yy()) +
	       minX.By());
  // Planer interpolation not required in z component
  Bfield[2] = minX.Bz();
  // The following is part of the attempted increase in efficiency. The
  //   idea was to carry the reference grid points forward to the next 
  //    iteration, by assigning them to variables in the class. The problem
  //   is that the assignment fails because of the const declaration of the 
  //   function 
  gridtested = true;
  _xmin = new double[6];
  _xmax = new double[6];
  _ymin = new double[6];
  _ymax = new double[6];
  _xmin[0] = minX.xx();_xmin[1] = minX.yy();_xmin[2] = minX.zz();
  _xmin[3] = minX.Bx();_xmin[4] = minX.By();_xmin[5] = minX.Bz();
  _xmax[0] = maxX.xx();_xmax[1] = maxX.yy();_xmax[2] = maxX.zz();
  _xmax[3] = maxX.Bx();_xmax[4] = maxX.By();_xmax[5] = maxX.Bz();
  _ymin[0] = minY.xx();_ymin[1] = minY.yy();_ymin[2] = minY.zz();
  _ymin[3] = minY.Bx();_ymin[4] = minY.By();_ymin[5] = minY.Bz();
  _ymax[0] = maxY.xx();_ymax[1] = maxY.yy();_ymax[2] = maxY.zz();
  _ymax[3] = maxY.Bx();_ymax[4] = maxY.By();_ymax[5] = maxY.Bz();
  
}
*/

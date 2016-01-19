#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <ostream>

#include <recpack/EVectorMap.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <mind/MINDfieldMapReader.h>
#include <recpack/Definitions.h>

//namespace RecPack {
using namespace std;
using namespace Recpack;

Recpack::MINDfieldMapReader::MINDfieldMapReader(const std::string& bf, double fieldscale)
  :EVectorMap(3), _fieldScale(fieldscale), _geom(1), _useMap(true)
{
  _Xmin =  999999.99;
  _Xmax = -999999.99;
  _Ymin =  999999.99;
  _Ymax = -999999.99;
  // _vecMap.reserve(m.size());
  std::map< double, std::map< double, EVector > > field_map;
  std::map<double, std::map<double, EVector> >::iterator ix;
  std::map<double, EVector>::iterator iy;
  ifstream bfile;
  bfile.open(bf.c_str(), ifstream::in);
  char temp[256];
  while( bfile.getline( temp, 256 ) ){
    std::string s(temp);
    std::vector<std::string> entries;
    Tokenize(s, entries);
    std::vector<double> points;
    for(std::vector<std::string>::iterator it = entries.begin();
	it != entries.end(); it++)
      points.push_back(atof((*it).c_str()));
    EVector B = EVector(3,0);
    B[0] = points[3]*tesla;
    B[1] = points[4]*tesla;
    B[2] = points[5]*tesla;
    // std::cout<<B[0]<<"\t"<<B[1]<<"\t"<<B[2]<<std::endl;
    field_map[points[0]*m][points[1]*m] = B;
  }
  bfile.close();
  // std::cout<<"Tesla scaling is "<<tesla<<std::endl;
  
  Nx = double(field_map.size());
  for( ix=field_map.begin(); ix != field_map.end(); ix++ ){
    if((*ix).first < _Xmin) _Xmin = (*ix).first;
    if((*ix).first > _Xmax) _Xmax = (*ix).first;
    std::vector<EVector> ytemp;
    if(ix == field_map.begin()) Ny = double((*ix).second.size());
    for( iy=(*ix).second.begin(); iy!=(*ix).second.end(); iy++){
      if((*iy).first < _Ymin) _Ymin = (*iy).first;
      if((*iy).first > _Ymax) _Ymax = (*iy).first;
      ytemp.push_back((*iy).second);
    }
    _vecMap.push_back(ytemp);
  }
  
  dx = (_Xmax - _Xmin)/Nx;
  dy = (_Ymax - _Ymin)/Ny;
  
}


Recpack::MINDfieldMapReader::MINDfieldMapReader(double fieldscale,
						int Geom, double detX, double detY)
  :EVectorMap(3), _fieldScale(fieldscale), _geom(fabs(Geom)), _useMap(false), 
   _Xmin(-detX/2.), _Xmax(detX/2.), _Ymin(-detY/2.), _Ymax(detY/2.)
{
  
}

EVector Recpack::MINDfieldMapReader::compute_vector(const EVector& pos) const {
  EVector BField = EVector(3,0);
  if(_useMap){
    if(pos[0] < _Xmin  || pos[0] > _Xmax || 
       pos[1] < _Ymin  || pos[1] > _Ymax){
      BField[0] = 0.0; // _vecMap[0][0][0];
      BField[1] = 0.0; // _vecMap[0][0][1];
      BField[2] = 0.0; // _vecMap[0][0][2];
      return BField;
    }
    else {
      // Assume a regular grid of points
      int j1 = int(floor((pos[0] - _Xmin)/dx + 0.5));
      int k1 = int(floor((pos[1] - _Ymin)/dy + 0.5));
      if(j1 < 0) j1 = 0;
      if(k1 < 0) k1 = 0;
      if(j1 >= _vecMap.size())     j1 = _vecMap.size() - 1;
      if(k1 >= _vecMap.at(j1).size()) k1 = _vecMap.at(j1).size() - 1;
      double x1 = _Xmin + double(j1)*dx;
      double x2 = _Xmin + double(j1-1)*dx;
      double y1 = _Ymin + double(k1)*dy;
      double y2 = _Ymin + double(k1-1)*dy;
      int j2 = j1-1;
      int k2 = k1-1;
      if ((pos[0] > x1 && j1 <= _vecMap.size() - 2)  || j1 == 0){
	x2 = _Xmin + double(j1 + 1)*dx;
	j2 = j1+1;
      }
      if ((pos[1] > y1 && k1 <= _vecMap[j1].size() - 2) || k1 == 0){
	y2 = _Ymin + double(k1 + 1)*dy;
	k2 = k1 + 1;
      }

      double Bx0 = _vecMap[j1][k1][0];
      double By0 = _vecMap[j1][k1][1];
      double Bx1 = _vecMap[j2][k1][0];
      double By1 = _vecMap[j2][k1][1];
      /*if( (Bx1 < 0.5*Bx0 || By1 < 0.5*By0) && j2 < j1 ){ 
	// A correction to avoid the edges of the plate
	j2 = j1+1;
	x2 = _Xmin + double(j2)*dx;
	Bx1 = _vecMap[j2][k1][0];
	By1 = _vecMap[j2][k1][1];
      }
      else if( (Bx1 < 0.5*Bx0 || By1 < 0.5*By0) && j2 > j1 ){
	j2 = j1-1;
	x2 = _Xmin + double(j2)*dx;
	Bx1 = _vecMap[j2][k1][0];
	By1 = _vecMap[j2][k1][1];
	}*/
      double Bx2 = _vecMap[j1][k2][0];
      double By2 = _vecMap[j1][k2][1];
      /*if( (Bx2 < 0.5*Bx0 || By2 < 0.5*By0) && k2 < k1 ){ 
	// A correction to avoid the edges of the plate
	k2 = k1+1;
	y2 = _Ymin + double(k2)*dx;
	Bx1 = _vecMap[j1][k2][0];
	By1 = _vecMap[j1][k2][1];
      }
      else if( (Bx2 < 0.5*Bx0 || By2 < 0.5*By0) && k2 > k1 ){
	k2 = k1-1;
	y2 = _Ymin + double(k2)*dx;
	Bx1 = _vecMap[j1][k2][0];
	By1 = _vecMap[j1][k2][1];
	}*/

      double dBxdx = (Bx1 - Bx0)/(x2 - x1);
      double dBydx = (By1 - By0)/(x2 - x1);
      double dBxdy = (Bx2 - Bx0)/(y2 - y1);
      double dBydy = (By2 - By0)/(y2 - y1);
      
      BField[0] = _fieldScale*(Bx0+dBxdx*(pos[0]-x1)+dBxdy*(pos[1]-y1));
      BField[1] = _fieldScale*(By0+dBydx*(pos[0]-x1)+dBydy*(pos[1]-y1));
      BField[2] = _fieldScale*(_vecMap[j1][k1][2]);
      
      return BField;
    }
  }
  else{
    if(_geom!=2){
      double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      //double th = atan2(pos[1],pos[0]);
      //double sinsq4th = pow(sin(4*th),2);
      BField[0] = -_fieldScale * ( 1.36 + 0.0406*m/r + 0.80*exp(-r*0.16/m) ) * pos[1]/r * tesla;
      BField[1] =  _fieldScale * ( 1.36 + 0.0406*m/r + 0.80*exp(-r*0.16/m) ) * pos[0]/r * tesla;
      BField[2] =  0.0 * tesla;
    
      return BField;
    } else {
      double rR = sqrt((pos[0] + _Xmax/2.)*(pos[0] + _Xmax/2.) + pos[1]*pos[1]);
      double rL = sqrt((pos[0] + _Xmin/2.)*(pos[0] + _Xmin/2.) + pos[1]*pos[1]);
      //double th = atan2(pos[1],pos[0]);
      //double sinsq4th = pow(sin(4*th),2);
      BField[0] = -_fieldScale * ( 1.36 + 0.0406*m/rL + 0.80*exp(-rL*0.16/m) ) * pos[1]/rL * tesla +
	_fieldScale * ( 1.36 + 0.0406*m/rR + 0.80*exp(-rR*0.16/m) ) * pos[1]/rR * tesla;
      BField[1] =  _fieldScale * ( 1.36 + 0.0406*m/rL + 0.80*exp(-rL*0.16/m) ) * (pos[0]+_Xmin/2.)/rL * tesla -
	_fieldScale * ( 1.36 + 0.0406*m/rR + 0.80*exp(-rR*0.16/m) ) * (pos[0]+_Xmax/2.)/rR * tesla;
      BField[2] =  0.0 * tesla;
      return BField;
    }
  }
}

  
void MINDfieldMapReader::Tokenize(const std::string& str,
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
//}

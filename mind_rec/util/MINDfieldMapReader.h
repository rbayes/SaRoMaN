#ifndef MINDfieldMapReader_h
#define MINDfieldMapReader_h 1

#include <recpack/EAlgebra.h>
#include <recpack/RecPackManager.h>
#include <recpack/EVectorMap.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>


using namespace Recpack;

namespace Recpack {
  //! A EVector depending on the position 
  class MINDfieldMapReader : public EVectorMap 
  {
  protected:
    
    //! the default vector
    EVector _vector;
    //! assume a longitudinally uniform field
    // std::map<double, std::map<double, EVector> > _VectorMap;
    //! generate a vector containing a rectangular grid of points from a map
    std::vector< std::vector< EVector > > _vecMap;
    double _fieldScale;
    int _geom;
    bool _useMap;
    double _Xmin;
    double _Xmax;
    double _Ymin;
    double _Ymax;
    double dx;
    double dy;
    double Nx;
    double Ny;
  private:
    void Tokenize(const std::string& str,
		  std::vector<std::string>& tokens,
		  const std::string& delimiters=" ");
    
  public:
    //! default constructor
    MINDfieldMapReader(size_t n=3){_vector = EVector(n,0);}
      
    //! constructer with a vector
    MINDfieldMapReader(EVector& v){_vector=v; }
  
    //! default destructor
    ~MINDfieldMapReader(){;}
  
    //! constructer with a file name in a const string
    MINDfieldMapReader(const std::string& bfile, double fieldscale);
    
    //! constructor with a field scaling
    MINDfieldMapReader(double fieldscale, int Geom, double detX, double detY);
    
    //! return non-const vector 
    EVector compute_vector(const EVector& pos) const;
    
    
  };
  
} 
#endif




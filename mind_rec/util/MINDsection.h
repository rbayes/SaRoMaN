#ifndef RPMINDsection_h
#define RPMINDsection_h 1

#include <recpack/Plane.h>

using namespace Recpack;

namespace Recpack{

  //! a cross-sectional face of the MIND detector
  class MINDsection : public Plane {	
  public:
    
    //! default constructor
    MINDsection(const EVector& pos,const EVector& normal,const EVector& uaxis,
		double U, double V, double fw, double fh);
    
    //! default destructor
    virtual ~MINDsection() {};
    

    //! return true if the point is inside
    bool is_inside(const EVector& ro, Messenger::Level level, 
		   double tolerance = EGeo::zero_l()) const;

  private :
    double _flange_width;
    double _flange_height;
    
  };
}
#endif

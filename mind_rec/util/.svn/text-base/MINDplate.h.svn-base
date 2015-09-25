#ifndef RPMINDplate_h
#define RPMINDplate_h 1

#include <recpack/Volume.h>
#include <recpack/Rectangle.h>
#include <MINDsection.h>

using namespace Recpack;

namespace Recpack{

  class MINDplate : public Volume {
  private:
    double fw;
    double fh;
  public:
    
    //! default constructor
    MINDplate(const EVector& pos, const EVector& uaxis, const EVector& vaxis, 
	      double U, double V, double W, 
	      double flange_width, double flange_height) 
      :Volume("MINDplate", pos)
      {create(uaxis,vaxis,U,V,W,flange_width,flange_height);}

    //! default destructor
    virtual ~MINDplate() {};

    //! return true if the point is inside
    bool is_inside(const EVector& r0) const;
    
  protected:
    //! create the plate of iron or scintillator
    void create(const EVector& uaxis, const EVector& vaxis,
		double U, double V, double W,
		double flange_width, double flange_height);
  };
}

#endif

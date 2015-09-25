#ifndef RPEMINDplate_h
#define RPEMINDplate_h 1

#include <recpack/Volume.h>
#include <recpack/Rectangle.h>
#include <MINDsection.h>

using namespace Recpack;

namespace Recpack{

  class EMINDplate : public Volume {
  private:
    double fw;
    double fh;
  public:
    
    //! default constructor
    EMINDplate(const EVector& pos, const EVector& uaxis, const EVector& vaxis, 
	      double U, double V, double W, 
	      double flange_width, double flange_height) 
      :Volume("EMINDplate", pos)
      {create(uaxis,vaxis,U,V,W,flange_width,flange_height);}

    //! default destructor
    virtual ~EMINDplate() {};

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

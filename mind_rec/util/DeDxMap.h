#ifndef DeDxMap_hxx
#define DeDxMap_hxx 1

#include <recpack/DynamicProperty.h>
#include <bhep/system_of_units.h>
//#include <HEPUnits.hxx>

using namespace Recpack;
using namespace bhep;
using namespace CLHEP;

#define NPBINS 28
const double  _de_dx_min_def = 0.195*MeV/mm;

namespace Recpack {

  //! A EVector depending on the position 
  class DeDxMap: public DynamicProperty {
  protected:

    double  _de_dx[NPBINS];
    double  _de_dx_elec[NPBINS];
    double  _pbins[NPBINS];
      
    double _melectron;
    double _mproton;
    double _mmuon;
    double _mpion;
    double _mkaon;

    double _de_dx_min;
  
  public:
    
    //! constructor with dimensions, null vector
    DeDxMap(double de_dx_min = _de_dx_min_def);
    
    //! default destructor
    virtual ~DeDxMap(){};

    //! return non-const vector 
    double property(const State& state) const;

    double de_dx_min(){return _de_dx_min;}
    
    // double getDeDx_map(double mom);///
  };
    
} 
#endif


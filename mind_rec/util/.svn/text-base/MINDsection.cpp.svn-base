#include <recpack/Rectangle.h>
#include <MINDsection.h>

//****************************************************
Recpack::MINDsection::MINDsection(const EVector& pos, const EVector& normal, const EVector& uaxis, double U, double V, double fw, double fh) : Plane(pos,normal) {
//****************************************************
  _names["shape"] = "mindsection";
  _axis["U"] = uaxis;
  _axis["V"] = crossprod(normal,uaxis);
  _parameters["U"] = fabs(U);
  _parameters["V"] = fabs(V);
  _flange_width  = fw;
  _flange_height = fh;
}

//****************************************************
bool Recpack::MINDsection::is_inside(const EVector& r0, Messenger::Level level,  
				   double tolerance) const {
//****************************************************

  EVector res = (r0-position());

  bool ok = true;
  const double Pi = 4*atan(1);
  double t45 = tan(Pi/8.);

  if (!Surface::is_inside(r0,level,tolerance)) ok = false;
  else if ( fabs(dot(res,_axis["V"])) > _parameters["V"] ) ok = false;
  else if ( fabs(dot(res,_axis["U"])) > 
	    _parameters["U"]* (1 + t45) - fabs(dot(res,_axis["V"])))
    ok = false;
  else if ( dot(res,_axis["V"]) > _flange_height ){
    if( fabs(dot(res,_axis["U"])) > _parameters["U"] + _flange_width)
      ok = false;
  }
  else if ( fabs(dot(res,_axis["U"])) > _parameters["U"] ) ok = false;


  if (Messenger::VERBOSE <= level){
    std::cout << "       Rectangle::is_inside(). ok = " << ok 
	      << " --> r = " << print(r0)  
	      << " *** l = (" << _parameters["U"] << ", "  << _parameters["V"] << ")" 
	      << " *** d = (" << dot(res,_axis["U"]) << ", " << dot(res,_axis["V"]) <<  ")" << std::endl;
  }

  //  if (Messenger::VVERBOSE <= level){
  //    std::cout << *this << std::endl;
  //  }

  return ok;
  
  
}


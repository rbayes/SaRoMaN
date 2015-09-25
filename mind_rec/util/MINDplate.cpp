#include <MINDplate.h>
#include <recpack/Rectangle.h>
#include <MINDsection.h>
#include <cmath>

//*******************************************
void Recpack::MINDplate::create(const EVector& uaxis, const EVector& vaxis,
				double U, double V, double W, 
				double flange_width, double flange_height) {
  //*****************************************


  const double Pi = 4*atan(1);
  double t22 = tan(Pi/8.);
  double s22 = 1./cos(Pi/8.);

  _axis["U"] = uaxis;
  _axis["V"] = vaxis;
  EVector waxis = crossprod(uaxis, vaxis);
  _axis["W"] = waxis;

  _parameters["U"] = fabs(U);
  _parameters["V"] = fabs(V);
  _parameters["W"] = fabs(W);
  _parameters["fw"] = flange_width;
  _parameters["fh"] = flange_height;

  double UU = fabs(U);
  double VV = fabs(V);
  double WW = fabs(W);
  
  fw = flange_width;
  fh = flange_height - UU;
  
  EVector deltaU = uaxis*UU;
  EVector deltaV = vaxis*VV;
  EVector deltaW = waxis*WW;
  double rt2 = sqrt(2.0);

  EVector deltaS = (uaxis * (UU + VV*t22)/2. + vaxis * (VV + UU*t22)/2.);
  EVector deltaT =(-uaxis * (UU + VV*t22)/2. + vaxis * (VV + UU*t22)/2.);
  double SS = 0.5*sqrt((UU*UU + VV*VV)*s22*s22 - 4*UU*VV*t22);

  EVector saxis = deltaS/sqrt(dot(deltaS,deltaS));
  EVector taxis = deltaT/sqrt(dot(deltaT,deltaT));

  EVector deltaQ = deltaS - VV/UU* fw/(rt2)*taxis;
  EVector deltaR = deltaT - VV/UU* fw/(rt2)*saxis;
  double QQ = sqrt(dot(deltaQ,deltaQ));
  double RR = sqrt(dot(deltaR,deltaR));
  EVector qaxis = 1./QQ * deltaQ;
  EVector raxis = 1./RR * deltaR;


  // S1 bottom
  add_surface("outer", new Rectangle(_position-deltaV,vaxis,waxis,WW,t22*VV));
  // S2 top
  add_surface("outer", new Rectangle(_position+deltaV,vaxis,waxis,WW,t22*VV));
  // S3 bottom-side -
  add_surface("outer", new Rectangle(_position-deltaS,saxis,waxis,WW,t22*UU));
  // S4 bottom-side +
  add_surface("outer", new Rectangle(_position-deltaT,taxis,waxis,WW,t22*UU));
  // S5 top-side +
  add_surface("outer", new Rectangle(_position+deltaQ,qaxis,taxis,
				     t22*UU + fw/rt2, WW));
  // S6 top-side -
  add_surface("outer", new Rectangle(_position+deltaR,raxis,saxis,
				     t22*UU + fw/rt2, WW));
  // S7 side - 
  add_surface("outer", new Rectangle(_position-deltaU-vaxis*(t22*UU-fh)/2.,
				     uaxis,vaxis,(t22*UU + fh)/2.0,WW));
  // S8 side + 
  add_surface("outer", new Rectangle(_position+deltaU-vaxis*(t22*UU-fh)/2.,
				     uaxis,vaxis,(t22*UU + fh)/2.0,WW));

  // S9 back
  add_surface("outer", new MINDsection(_position-deltaW,waxis,uaxis,UU,VV,fw,fh));
  // S10 front
  add_surface("outer", new MINDsection(_position+deltaW,waxis,uaxis,UU,VV,fw,fh));
  
  // S11 flange tip +
  add_surface("outer", new Rectangle(_position+ uaxis*(UU + fw) + 
				     vaxis*0.5*(UU*t22 - fw*VV/UU + fh), 
				     uaxis, vaxis,0.5*(UU*t22 - fw*VV/UU - fh),
				     WW));
  // S12 flange tip -
  add_surface("outer", new Rectangle(_position- uaxis*(UU + fw) + 
				     vaxis*0.5*(UU*t22 - fw*VV/UU + fh), 
				     uaxis, vaxis,0.5*(UU*t22 - fw*VV/UU - fh),
				     WW));
  // S13 flange bottom +
  add_surface("outer", new Rectangle(_position + uaxis*(UU + fw/2.) + 
				     vaxis*fh, vaxis, uaxis, fw/2., WW));
  // S13 flange bottom -
  add_surface("outer", new Rectangle(_position - uaxis*(UU + fw/2.) + 
				     vaxis*fh, vaxis, uaxis, fw/2., WW));
  // _size = uaxis*(UU + fw) + vaxis*VV + waxis*WW;

}

bool Recpack::MINDplate::is_inside(const EVector& r0) const {
  
  EVector res = (r0-position());
  
  const double Pi = 4*atan(1);
  double t22 = tan(Pi/8.);
  double slope = (_parameters["V"] - t22*_parameters["U"])/
    (_parameters["U"]*t22 - _parameters["V"]);

  if ( fabs(dot(res,_axis["V"])) >
       _parameters["U"]*t22 + slope*(fabs(dot(res,_axis["U"])) - 
				     _parameters["U"]))
       return false;
  if ( fabs(dot(res,_axis["V"])) > _parameters["V"] ) return false;
  if ( dot(res,_axis["V"]) > _parameters["fh"] + _parameters["V"] ){
    if( fabs(dot(res,_axis["U"])) > _parameters["U"] + _parameters["fw"])
      return  false;
  else if ( fabs(dot(res,_axis["U"])) > _parameters["U"] ) return false; }
  if (fabs(dot(res,_axis["W"])) > _parameters["W"]) return false;
  
  return true;
  
}

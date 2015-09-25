#include <rec_hit.h>

rec_hit::rec_hit(): Measurement()
{
  _totEng = 0;
  _muProp = 0;
}

rec_hit::~rec_hit()
{
}

void rec_hit::add_point(bhep::hit *point)
{
  //Add the point to the contributing vector.

  _points.push_back( point );

  _totEng += point->ddata("EnergyDep");

  string moth_name = point->sdata( "true_moth" );
  int moth_prim = point->idata( "moth_prim" );

  if ( moth_name == "mu+" || moth_name == "mu-" ){
    if ( moth_prim == 1 )
      _muProp++;
  } else if ( moth_name == "lepton_shower" )
    _muProp += 0.5;

  _npoints = _points.size();

}

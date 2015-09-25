#include <shower_rec.h>

#include <bhep/utilities.h>

//*******************************
shower_rec::shower_rec(){
//*******************************

}

//*******************************
bool shower_rec::initialize( double tol, double resolution ){
//*******************************

  _tolerance = tol * cm;
  _res = resolution * cm / sqrt( 12 );//????

  _dirVec = EVector(3,0);

  return true;
}

//*******************************
bool shower_rec::execute( shower& had_shower ){
//*******************************
  cout<<"execute shower_rec "<<endl;
  bool ok;

  reset();

  order_hits( had_shower.get_had_hits() );
 
  ok = reconstruct( had_shower.get_vertex() );

  if ( ok )
    had_shower.set_unit_direction( _dirVec );

  return ok;
}

//*******************************
bool shower_rec::finalize(){
//*******************************

  return true;
}

//*******************************
void shower_rec::reset(){
//*******************************

  _zpos.clear();
  for (unsigned int ii = 0;ii < _planes;ii++)
    _meas_in_plane[ii].clear();
  _planes = 0;
  _dirVec[0] = 0; _dirVec[1] = 0; _dirVec[2] = 0;
  
}
//*******************************
void shower_rec::order_hits( measurement_vector& hads ){
//*******************************
//Start by putting the hits in order or ascending Z.
  
  sort( hads.begin(), hads.end(), forwardSorter() );

  measurement_vector::iterator hadIt;
  size_t nplanes = 0;
  double z_current;

  _meas_in_plane[nplanes].push_back( ( *(hads.begin()) ) );
  _zpos.push_back( ( *(hads.begin()) )->position()[2] );
  double z_prev = _zpos[nplanes];
  
  for (hadIt = hads.begin()+1;hadIt != hads.end();hadIt++){
    
    z_current = (*hadIt)->position()[2];
    
    if ( fabs( z_current - z_prev ) >= _tolerance ){
      nplanes++;
      _zpos.push_back( z_current );
    } else _zpos[nplanes] += (*hadIt)->position()[2];
    
    _meas_in_plane[nplanes].push_back( (*hadIt) );
    
    z_prev = z_current;
    
  }
  
  for (unsigned int inorm = 0;inorm <= nplanes;inorm++)
    _zpos[inorm] /= (double)_meas_in_plane[inorm].size();
  
  _planes = nplanes+1;
  
}

//*******************************
bool shower_rec::reconstruct( EVector vert ){
//*******************************
//currently a very simple straight line fit
//which does not consider vertex.

  if ( _planes < 2 ) return false;

  measurement_vector::iterator plIt;
  const dict::Key Edep = "E_dep";

  double x[_planes], y[_planes], z[_planes];
  double xerr[_planes], yerr[_planes], zerr[_planes];
  double hit_eng, eng_plane;

  for (int iplane = 0;iplane < (int)_planes;iplane++){

    x[iplane] = 0;
    y[iplane] = 0;
    z[iplane] = _zpos[iplane];
    xerr[iplane] = _res;
    yerr[iplane] = _res;
    zerr[iplane] = 0.;
    eng_plane = 0;
    
    for (plIt = _meas_in_plane[iplane].begin();plIt != _meas_in_plane[iplane].end();plIt++)
      {
	hit_eng = bhep::double_from_string( (*plIt)->name( Edep ) ) * GeV;
	
	eng_plane += hit_eng;
	
	x[iplane] += (*plIt)->position()[0] * hit_eng;
	y[iplane] += (*plIt)->position()[1] * hit_eng;

      }
    
    x[iplane] /= eng_plane;
    y[iplane] /= eng_plane;
    
  }
  
  TGraphErrors * gr1 = new TGraphErrors( _planes, z, x, zerr, xerr );
  TGraphErrors * gr2 = new TGraphErrors( _planes, z, y, zerr, xerr );

  TF1 *fitfunc = new TF1("fitfun","[0]+[1]*x",-3,3);

  gr1->Fit("fitfun","QN");
  _dirVec[0] = fitfunc->GetParameter(1);

  gr2->Fit("fitfun","QN");
  _dirVec[1] = fitfunc->GetParameter(1);

  _dirVec[2] = 1;

  _dirVec /= _dirVec.norm();

  delete gr1;
  delete gr2;
  delete fitfunc;

  return true;
}

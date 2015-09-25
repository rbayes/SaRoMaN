#include <hit_constructor.h>

// #include <CLHEP/Random/RandGauss.h>
#include <TRandom3.h>

hit_constructor::hit_constructor(const bhep::gstore& store)
{
  //Store detector geometry and calculate scint layer positions.

  _detectorLength = store.fetch_dstore("MIND_z") * m;
  _detectorX = store.fetch_dstore("MIND_x") * m;
  _detectorY = store.fetch_dstore("MIND_y") * m;
  _passiveLength = store.fetch_dstore("widthI") * cm;
  _activeLength = store.fetch_dstore("widthS") * cm;
  _braceLength = 0.0;
  if (store.find_dstore("widthAl"))
    _braceLength = store.fetch_dstore("widthAl") * cm;
  _gapLength = store.fetch_dstore("widthA") * cm;
  _nActive = store.fetch_istore("nplane");
  OctGeom = 0;
  if (store.find_istore("isOctagonal"))
    OctGeom = store.fetch_istore("isOctagonal");
  _voxXdim = store.fetch_dstore("rec_boxX") * cm;
  _voxYdim = store.fetch_dstore("rec_boxY") * cm;
  _nVoxX = (int)( _detectorX / _voxXdim );
  _nVox = _nVoxX * (int)( _detectorY / _voxYdim );

  long seed = (long)store.fetch_dstore("Gen_seed");
  _ranGen = TRandom3( seed );

  _minEng = store.fetch_dstore("min_eng") * MeV;

  // double res = store.fetch_dstore("pos_res") * cm;
//   _cov = EMatrix(2,2,0);
//   for (int i = 0;i < 2;i++)
//     _cov[i][i] = pow(res, 2);

  //_measType = store.fetch_sstore("meas_type");

  calculate_layerZ();
  
}

hit_constructor::~hit_constructor()
{
}

void hit_constructor::reset()
{
  //Clear out map.correct??
  _voxels.clear();

}

void hit_constructor::execute(const std::vector<bhep::hit*>& hits,
			      std::vector<bhep::hit*>& rec_hit)
{
  //First clear out map.
  reset();
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;

  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );
  
  //sort into voxels map.
  parse_to_map( sortedHits );
  
  //Make rec_hits from vox.
  construct_hits( rec_hit );
  
}

void hit_constructor::calculate_layerZ()
{
  //Fill vector with all possible scint z positions.

  double pieceLength = _passiveLength + _nActive*_activeLength
    + (_nActive+1)*_gapLength + (_nActive * _braceLength );
  int npieces = (int)ceil( _detectorLength / pieceLength );

  //reset detector length to integer multiple of pieces.
  _detectorLength = npieces * pieceLength;

  int nLayers = npieces * _nActive;

  double z;

  for (int iLayer = 0;iLayer < nLayers;iLayer++){

    z = 0;

    for (int j = 0;j < _nActive;j++){

      z += _gapLength + _activeLength;

      if ( (iLayer-j) % _nActive == 0 ){
	int block_no = (iLayer-j) / _nActive;
	z += pieceLength*block_no + _passiveLength
	  - _activeLength/2 + _braceLength;
	break;
      }

    }

    z -= _detectorLength/2;
    
    _zLayer.push_back( z );

  }
  
}

void hit_constructor::parse_to_map(const std::vector<bhep::hit*> hits)
{
  //Sort hits into voxel map.
  int voxNo;
  double zpos;

  //Set iterator to plane 1.
  _zIt = _zLayer.begin();

  std::vector<bhep::hit*>::const_iterator hitIt;

  for (hitIt = hits.begin();hitIt != hits.end();hitIt++){

    //find plane.
    zpos = find_plane( *(*hitIt) );
    
    //get Vox number;
    voxNo = calculate_vox_no( *(*hitIt) );

    //Add voxel to map (or hit to existing voxel).

    _voxels[zpos].insert( pair<int,bhep::hit*>(voxNo,(*hitIt)) );

  }
  
}

double hit_constructor::find_plane(bhep::hit& curHit)
{
  //find the appropriate z position by comparison to
  //Layer z.
  
  double modDiff = fabs( curHit.x()[2] - (*_zIt) );
  
  while ( (int)modDiff > (int)_activeLength/2 ){

    _zIt++;

    modDiff = fabs( curHit.x()[2] - (*_zIt) );
    
  }
  
  return (*_zIt);
}

int hit_constructor::calculate_vox_no(bhep::hit& curHit)
{
  //calculate the correct voxel number.

  int xbox = (int)( fabs( curHit.x()[0] + _detectorX/2 ) / _voxXdim );

  int ybox = (int)( fabs( curHit.x()[1] - _detectorY/2 ) / _voxYdim );

  int vox_num = xbox + ybox*_nVoxX;

  return vox_num;
}

void hit_constructor::construct_hits(std::vector<bhep::hit*>& rec_hit)
{
  //takes the voxels which have been filled and make
  //rec_hit objects out of them.

  std::map<double,std::multimap<int,bhep::hit*> >::iterator vIt;
  std::multimap<int,bhep::hit*>::iterator vIt2;

  for (vIt = _voxels.begin();vIt != _voxels.end();vIt++){

    while ( vIt->second.size() != 0 ){
      //set second iterator to first filled voxel in layer.
      vIt2 = vIt->second.begin();

      bhep::hit* vhit = get_vhit( vIt2->first, vIt->first, vIt->second );
      
      if ( vhit != NULL )
	rec_hit.push_back( vhit );

      vIt->second.erase( vIt2->first );

    }
  }
  
}

bhep::hit* hit_constructor::get_vhit(int vox, double z,
				     const std::multimap<int,bhep::hit*>& map1)
{
  //Makes a rec_hit from the voxel position and adds the relevant points.
  
  double totEng = 0., muProp = 0.; //these will be done in rec_hit eventually.
  vdouble X, Y, Z, E, T; //annoying but again all in rec_hit class.
  double proptime=9999999.9, vlight = 299792458. / 1.6;
  double meanvoxtime;
  int irow = vox / _nVoxX;
  int icol = vox % _nVoxX;
  
  double voxX = icol*_voxXdim + _voxXdim/2 - _detectorX/2;
  double voxY = _detectorY/2 - (irow*_voxYdim + _voxYdim/2);
  

  Point3D hitPos( voxX, voxY, z );
  
  bhep::hit* vhit = new bhep::hit( "tracking" );
  vhit->set_point( hitPos );

  vhit->add_property( "voxel", vox );

  std::multimap<int,bhep::hit*>::const_iterator hIt;
  for (hIt = map1.equal_range(vox).first;hIt != map1.equal_range(vox).second;hIt++)
    {
      X.push_back( (*hIt).second->x()[0] );
      Y.push_back( (*hIt).second->x()[1] );
      Z.push_back( (*hIt).second->x()[2] );
      T.push_back( (*hIt).second->ddata( "time" ) );
      E.push_back( (*hIt).second->ddata( "EnergyDep" ) );
      totEng += (*hIt).second->ddata( "EnergyDep" );
      proptime = (*hIt).second->ddata( "time" )  < proptime ?  
	(*hIt).second->ddata( "time" )  : proptime;
      if ( (*hIt).second->mother_particle().name() == "mu+" ||
	   (*hIt).second->mother_particle().name() == "mu-" ){
	if ( (*hIt).second->mother_particle().fetch_sproperty("CreatorProcess")=="none" )
	  muProp++;
      } else if ( (*hIt).second->mother_particle().name() == "lepton_shower" )
	muProp += 0.5;
      
      //std::cout<<(*hIt).second->x()[0]<<"\t"<<(*hIt).second->x()[1]<<"\t"<<(*hIt).second->x()[2]<<"\t"
      // <<(*hIt).second->ddata( "EnergyDep" )<<std::endl;
    }
  //Do attenuations.
  double xE1, xE2, yE1, yE2;
  double Xphot, Yphot;
  //Assume equal amounts of energy from both views and
  // equal energy flow in both directions along strip.
  xE1 = xE2 = yE1 = yE2 = totEng/4;
  //   proptime /= T.size();
  //   std::cout<<proptime<<std::endl;

  double slope = OctGeom == 1 ? (_detectorY - _detectorX*tan(atan(1)/2.))/
    (_detectorY*tan(atan(1)/2.) - _detectorX) : -1.;
  double dt, dtx, dty;
  //need to take into account drift distance to closest and furthest edge.
  if(!OctGeom){ // Assume a cylindrical geometry
    double xedge = _detectorX * sqrt(1 - pow(voxY/_detectorY, 2));
    double yedge = _detectorY * sqrt(1 - pow(voxX/_detectorX, 2));
    xE1 = xE1 * exp( -(xedge - fabs(voxX))/5000 );
    xE2 = xE2 * exp( -(xedge + fabs(voxX))/5000 );
    yE1 = yE1 * exp( -(yedge - fabs(voxY))/5000 );
    yE2 = yE2 * exp( -(yedge + fabs(voxY))/5000 );
    dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
      (xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
    dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
      (yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
    dt = dtx < dty ? dtx : dty;
  }
  else if(OctGeom==1){
    if((fabs(voxX) < _detectorY*tan(atan(1)/2.)/2 &&
	fabs(voxY) < _detectorX*tan(atan(1)/2.)/2 ) ){
      xE1 = xE1 * exp( -(_detectorX/2 - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(_detectorX/2 + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(_detectorY/2 - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(_detectorY/2 + fabs(voxY))/5000 );
      double xedge = _detectorX/2., yedge = _detectorY/2.;
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else if(fabs(voxX) > _detectorY*tan(atan(1)/2.)/2 &&
	    fabs(voxY) < _detectorX*tan(atan(1)/2.)/2  ){
      double xedge = _detectorX/2 + (fabs(voxY) - 
				     _detectorX/2.*tan(atan(1)/2.))*slope;
      xE1 = xE1 * exp( -(xedge  - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(xedge + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(_detectorY/2 - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(_detectorY/2 + fabs(voxY))/5000 );
      double yedge = _detectorY/2.;
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else if(fabs(voxX) < _detectorY*tan(atan(1)/2.)/2 &&
	    fabs(voxY) > _detectorX*tan(atan(1)/2.)/2  ){
      double yedge = _detectorY/2 + (fabs(voxX) - 
				     _detectorY/2.*tan(atan(1)/2.))*slope;
      xE1 = xE1 * exp( -(_detectorX/2 - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(_detectorX/2 + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(yedge - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(yedge + fabs(voxY))/5000 );
      double xedge = _detectorX/2.;
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else if(fabs(voxX) > _detectorY*tan(atan(1)/2.)/2 &&
	    fabs(voxY) > _detectorX*tan(atan(1)/2.)/2  ){
      double xedge = _detectorX/2 + (fabs(voxY) - 
				     _detectorX/2.*tan(atan(1)/2.))*slope;
      double yedge = _detectorY/2 + (fabs(voxX) - 
				     _detectorY/2.*tan(atan(1)/2.))*slope;
      xE1 = xE1 * exp( -(xedge - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(xedge + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(yedge - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(yedge + fabs(voxY))/5000 );
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
  }else if(OctGeom==2){
    
    if((fabs(voxX) < _detectorY*tan(atan(1)/2.)/2 + _detectorX/4. &&
	fabs(voxY) < _detectorY*tan(atan(1)/2.)/2 ) ){
      xE1 = xE1 * exp( -(_detectorX/2 - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(_detectorX/2 + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(_detectorY/2 - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(_detectorY/2 + fabs(voxY))/5000 );
      double xedge = _detectorX/2., yedge = _detectorY/2.;
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else if(fabs(voxX) > _detectorY*tan(atan(1)/2.)/2 + _detectorX/4. &&
	    fabs(voxY) < _detectorY*tan(atan(1)/2.)/2  ){
      double xedge = _detectorX/2 + (fabs(voxY) - 
				     _detectorY/2.*tan(atan(1)/2.))*slope;
      xE1 = xE1 * exp( -(xedge  - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(xedge + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(_detectorY/2 - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(_detectorY/2 + fabs(voxY))/5000 );
      double yedge = _detectorY/2.;
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else if(fabs(voxX) < _detectorY*tan(atan(1)/2.)/2 + _detectorX/4.&&
	    fabs(voxY) > _detectorY*tan(atan(1)/2.)/2  ){
      double yedge = _detectorX/2 + (fabs(voxX) - 
				     _detectorY/2.*tan(atan(1)/2.))*slope;
      xE1 = xE1 * exp( -(_detectorX/2 - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(_detectorX/2 + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(yedge - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(yedge + fabs(voxY))/5000 );
      double xedge = _detectorX/2.;
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else if(fabs(voxX) > _detectorY*tan(atan(1)/2.)/2  + _detectorX/4.&&
	    fabs(voxY) > _detectorX*tan(atan(1)/2.)/2  ){
      double xedge = _detectorX/2 + (fabs(voxY) - 
				     _detectorY/2.*tan(atan(1)/2.))*slope;
      double yedge = _detectorX/2 + (fabs(voxX) - 
				     _detectorY/2.*tan(atan(1)/2.))*slope;
      xE1 = xE1 * exp( -(xedge - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(xedge + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(yedge - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(yedge + fabs(voxY))/5000 );
      
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
  }
  else if(OctGeom == -2){ // Assume a square cross-section
    double xedge = _detectorY/2.;
    double yedge = _detectorY/2.;
    if((fabs(voxX) < _detectorX/2. &&
	fabs(voxY) < _detectorY/2. ) ){
      xE1 = xE1 * exp( -(xedge - fabs(voxX))/5000 );
      xE2 = xE2 * exp( -(xedge + fabs(voxX))/5000 );
      yE1 = yE1 * exp( -(yedge - fabs(voxY))/5000 );
      yE2 = yE2 * exp( -(yedge + fabs(voxY))/5000 );
      dtx = (xedge - fabs(voxX)) < (xedge + fabs(voxX)) ?
	(xedge - fabs(voxX))/vlight : (xedge + fabs(voxX))/vlight;
      dty = (yedge - fabs(voxY)) < (yedge + fabs(voxY)) ?
	(yedge - fabs(voxY))/vlight : (yedge + fabs(voxY))/vlight;
      dt = dtx < dty ? dtx : dty;
    }
    else {
      xE1 = 0.;
      xE2 = 0.;
      yE1 = 0.;
      yE2 = 0.;
    }
  }

  //smear the reconstructed energies.
  double sigMaX1 = 0.06 * xE1;
  double sigMaX2 = 0.06 * xE2;
  double sigMaY1 = 0.06 * yE1;
  double sigMaY2 = 0.06 * yE2;
  //smear and recombine.
  double xE = xE1 + _ranGen.Gaus( 0, sigMaX1 )
    + xE2 + _ranGen.Gaus( 0, sigMaX2 );
  double yE = yE1 + _ranGen.Gaus( 0, sigMaY1 )
    + yE2 + _ranGen.Gaus( 0, sigMaY2 );
  
  if ( xE < _minEng && yE < _minEng ){
    delete vhit;
    return NULL;
  }

  totEng = xE + yE;
  proptime += dt;

  vhit->add_property( "TotalEng", totEng );
  vhit->add_property( "XEng", xE );
  vhit->add_property( "YEng", yE );
  vhit->add_property( "MuonProportion", muProp );
  vhit->add_property( "NoPoints", (int)X.size() );
  vhit->add_property( "Xpoint", X );
  vhit->add_property( "Ypoint", Y );
  vhit->add_property( "Zpoint", Z );
  vhit->add_property( "Epoint", E );
  vhit->add_property( "HitTime", proptime);

  return vhit;
}

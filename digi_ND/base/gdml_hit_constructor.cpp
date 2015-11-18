#include <gdml_hit_constructor.h>
#include <limits>

// #include <CLHEP/Random/RandGauss.h>
#include <TRandom3.h>
#include <algorithm>
gdml_hit_constructor::gdml_hit_constructor(const bhep::gstore& store)
{
  //Store detector geometry and calculate scint layer positions.

  _detectorLength = store.fetch_dstore("MIND_z") * m;
  _detectorX = store.fetch_dstore("MIND_x") * m;
  _detectorY = store.fetch_dstore("MIND_y") * m;
  _vertexDetdepth = store.fetch_dstore("vertex_z") * m;
  _vertexDetX = store.fetch_dstore("vertex_x") * m;
  _vertexDetY = store.fetch_dstore("vertex_y") * m;
  _passiveLength = store.fetch_dstore("passive_thickness") * cm;
  _activeLength = store.fetch_dstore("active_thickness") * cm;
  _braceLength = 0.0;
  if (store.find_dstore("bracing_thickness"))
    _braceLength = store.fetch_dstore("bracing_thickness") * cm;
  _gapLength = store.fetch_dstore("air_gap") * cm;
  _nActive = store.fetch_istore("active_layers");
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

  _attLength = store.fetch_dstore("WLSatten");

  calculate_layerZ();
  
}

gdml_hit_constructor::~gdml_hit_constructor()
{
}
/*
void gdml_hit_constructor::reset()
{
  //Clear out map.correct??
  _voxels.clear();

}
*/
void gdml_hit_constructor::execute(const std::vector<bhep::hit*>& hits,
			      std::vector<bhep::hit*>& rec_hit)
{
  //First clear out map.
  reset();
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;

  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );
  
  //sort into voxels map.
  //parse_to_map( sortedHits );
  
  //Make rec_hits from vox.
  //construct_hits( rec_hit );

  construct_hits(sortedHits, rec_hit);
  
}

void gdml_hit_constructor::calculate_layerZ()
{
  // should read this from root or gdml.

  
  //Fill vector with all possible scint z positions.
  double vertModlength = _nActive*_activeLength + _nActive*_gapLength;
  int nVertPieces = (int)ceil(_vertexDetdepth/vertModlength);
  _vertexDetdepth = nVertPieces * vertModlength;
  
  int nVertLayers = nVertPieces *_nActive;

  double pieceLength = _passiveLength + _nActive*_activeLength
    + (_nActive)*_gapLength + (_nActive * _braceLength );
  int npieces = (int)ceil( _detectorLength / pieceLength );

  //reset detector length to integer multiple of pieces.
  _detectorLength = npieces * pieceLength;

  int nLayers = npieces * _nActive;

  double z;
  
  for (int iLayer = 0;iLayer < nVertLayers;iLayer++){
    
    z = 0;
    
    for (int j = 0;j < _nActive;j++){

      z += _gapLength + _activeLength;

      if ( (iLayer-j) % _nActive == 0 ){
	int block_no = (iLayer-j) / _nActive;
	z += vertModlength*block_no 
	  - _activeLength/2;
	break;
      }

    }
    z -= (_detectorLength + _vertexDetdepth)/2.;

    _zLayer.push_back( z );
  }

  for (int iLayer = 0;iLayer < nLayers;iLayer++){

    z = 0;

    for (int j = 0;j < _nActive;j++){

      z += _gapLength + _activeLength;

      if ( (iLayer-j) % _nActive == 0 ){
	int block_no = (iLayer-j) / _nActive;
	z += pieceLength*block_no + _passiveLength
	  - _activeLength/2;
	break;
      }

    }

    z -= _detectorLength/2;
    z += _vertexDetdepth/2;
    
    _zLayer.push_back( z );

  }
  
}
/*
void gdml_hit_constructor::parse_to_map(const std::vector<bhep::hit*> hits)
{
  //Sort hits into voxel map.
  int voxNo;
  double zpos;

  //Set iterator to plane 1.
  //_zIt = _zLayer.begin();

  std::vector<bhep::hit*>::const_iterator hitIt;

  for (hitIt = hits.begin();hitIt != hits.end();hitIt++){
    
    _zIt = _zLayer.begin();
    //find plane.
    zpos = find_plane( *(*hitIt) );
    
    //get Vox number;
    voxNo = calculate_vox_no( *(*hitIt) );

    //Add voxel to map (or hit to existing voxel).
    if ( voxNo >= 0)
      _voxels[zpos].insert( pair<int,bhep::hit*>(voxNo,(*hitIt)) );
    
  }
  
}
*/
double gdml_hit_constructor::find_plane(bhep::hit& curHit)
{
  //find the appropriate z position by comparison to
  //Layer z.
  
  double modDiff = fabs( curHit.x()[2] - (*_zIt) );
  std::cout<<"X = "<<curHit.x()[0]<<", Y = "<<curHit.x()[1]<<", Z of hit "<<curHit.x()[2]<<std::endl;
  while ( (int)modDiff > (int)_activeLength/2. && _zIt != _zLayer.end()){

    _zIt++;

    modDiff = fabs( curHit.x()[2] - (*_zIt) );
    
  }
  
  return (*_zIt);
}
/*
int gdml_hit_constructor::calculate_vox_no(bhep::hit& curHit)
{
  //calculate the correct voxel number.
  int vox_num = -1;
  if ( fabs(curHit.x()[0]) < _detectorX/2 && fabs(curHit.x()[1]) < _detectorY/2 ){
    int xbox = (int)( fabs( curHit.x()[0] + _detectorX/2 ) / _voxXdim );
    
    int ybox = (int)( fabs( curHit.x()[1] - _detectorY/2 ) / _voxYdim );
    
    vox_num = xbox + ybox*_nVoxX;
  }    
	 
  return vox_num;
}
*/

void gdml_hit_constructor::construct_hits(const std::vector<bhep::hit*>& hits std::vector<bhep::hit*>& rec_hit)
{
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;
  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );

  // For each (sorted) hit, take the x,y,z positions smear these given the smearing and attenuation.

  for (std::vector<bhep::hit*>::iterator sortedHitIter = sortedHits.begin() ; sortedHitIter != sortedHits.end(); ++sortedHitIter)
    {
      bhep::hit* vhit = get_vhit(sortedHitIter->x()[0],sortedHitIter->x()[1],sortedHitIter->x()[2]);

      if ( vhit != NULL )
	{
	  rec_hit.push_back( vhit );
	}

    }



}

 /*
void gdml_hit_constructor::construct_hits(std::vector<bhep::hit*>& rec_hit)
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
      
      if ( vhit != NULL ){
	rec_hit.push_back( vhit );
	
	std::cout<<"x = "<<vhit->x()[0]<<",y = "<<vhit->x()[1]<<", z = "<<vhit->x()[2]<<std::endl;
      }
      vIt->second.erase( vIt2->first );
    }
  }
  
}

 */

bhep::hit* gdml_hit_constructor::get_vhit(double voxX,double voxY, double z,
				     const std::multimap<int,bhep::hit*>& map1)
{
  //Makes a rec_hit from the voxel position and adds the relevant points.
  bhep::hit* returnPointer;
  
  double totEng = 0., muProp = 0.; //these will be done in rec_hit eventually.
  vdouble X, Y, Z, E, T; //annoying but again all in rec_hit class.
  double proptime=9999999.9, vlight = 299792458. / 1.6;
  //double proptime = numeric_limits<double>::max()
  double meanvoxtime;
  int irow = vox / _nVoxX;
  int icol = vox % _nVoxX;

  double xedge = _detectorX/2.;
  double yedge = _detectorY/2.;
  double smearingFactor = 0.06;
  
  //double voxX = icol*_voxXdim + _voxXdim/2 - _detectorX/2;
  //double voxY = _detectorY/2 - (irow*_voxYdim + _voxYdim/2);
  
  Point3D hitPos( voxX, voxY, z );
  
  bhep::hit* vhit = new bhep::hit( "tracking" );
  vhit->set_point( hitPos );

  vhit->add_property( "voxel", vox );
  /*
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
      //	       <<(*hIt).second->ddata( "EnergyDep" )<<std::endl;
    }
*/

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

  xE1 = xE1 * exp(-(xedge - voxX)/_attLength);
  xE2 = xE2 * exp(-(3*xedge-voxX)/_attLength);
  yE1 = yE1 * exp(-(yedge - voxY)/_attLength);
  yE2 = yE2 * exp(-(3*yedge-voxY)/_attLength);
  
  dtx = (xedge - voxX) < (3*xedge - voxX) ?
    (xedge - voxX)/vlight : (3*xedge - voxX)/vlight;
  dty = (yedge - voxY) < (3*yedge - voxY) ?
    (yedge - voxY)/vlight : (3*yedge - voxY)/vlight;

  dt = dtx < dty ? dtx : dty;
  
  //smear the reconstructed energies.
  //smear and recombine.
  
  double xE = xE1 + _ranGen.Gaus( 0, smearingFactor * xE1 )
    + xE2 + _ranGen.Gaus( 0, smearingFactor * xE2 );
  double yE = yE1 + _ranGen.Gaus( 0, smearingFactor * yE1 )
    + yE2 + _ranGen.Gaus( 0, smearingFactor * yE2 );
  
  if ( fabs(z) > (_detectorLength + _vertexDetdepth)/2. && 
       xE < _minEng && yE < _minEng )
    {
      delete vhit;
      returnPointer = NULL;
    }
  else
    {
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
      
      returnPointer = vhit;
    }
  
  return returnPointer;
}

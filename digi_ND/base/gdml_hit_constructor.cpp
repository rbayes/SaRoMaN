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

  // calculate_layerZ();
  
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
  //reset();

  cout<<"In execute in gdml_hit_constructor.cpp"<<endl;
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;

  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );
  
  //sort into voxels map.
  //parse_to_map( sortedHits );
  
  //Make rec_hits from vox.
  //construct_hits( rec_hit );

  construct_hits(sortedHits, rec_hit);
  
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

void gdml_hit_constructor::construct_hits(const std::vector<bhep::hit*>& hits, std::vector<bhep::hit*>& rec_hit)
{

cout<<"In construct_hits in gdml_hit_constructor.cpp"<<endl;

  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;
  cout<<hits.size()<<endl;
  cout<<sortedHits.size()<<endl;
  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );

  // For each (sorted) hit, take the x,y,z positions smear these given the smearing and attenuation.

  for (std::vector<bhep::hit*>::iterator sortedHitIter = sortedHits.begin() ; sortedHitIter != sortedHits.end(); ++sortedHitIter)
    {

      cout<<"In for loop in execute in gdml_hit_constructor.cpp"<<endl;

      bhep::hit* vhit = get_vhit((*sortedHitIter));

      if ( vhit != NULL )
	{
	  rec_hit.push_back( vhit );
	}

    }
}


 bhep::hit* gdml_hit_constructor::get_vhit(bhep::hit* curr_hit)
{
cout<<"In get_vhit in gdml_hit_constructor.cpp"<<endl;

  double x = curr_hit->x()[0];
  double y = curr_hit->x()[1];
  double z = curr_hit->x()[2];


  //Makes a rec_hit from the voxel position and adds the relevant points.
  bhep::hit* returnPointer;
  
  double totEng = 0., muProp = 0.; //these will be done in rec_hit eventually.
  vdouble X, Y, Z, E, T; //annoying but again all in rec_hit class.
  double proptime=9999999.9, vlight = 299792458. / 1.6;
  //double proptime = numeric_limits<double>::max()
  //double meanvoxtime;
  //int irow = vox / _nVoxX;
  //int icol = vox % _nVoxX;

  double xedge = _detectorX/2.;
  double yedge = _detectorY/2.;
  double smearingFactor = 0.06;
  
  Point3D hitPos( x, y, z );
  
  bhep::hit* vhit = new bhep::hit( "tracking" );
  vhit->set_point( hitPos );

  X.push_back(x);
  Y.push_back(y);
  Z.push_back(z);
  T.push_back(curr_hit->ddata( "time" ) );
  E.push_back(curr_hit->ddata( "EnergyDep" ) );
  totEng += curr_hit->ddata( "EnergyDep" );
  proptime = curr_hit->ddata( "time" )  < proptime ?  
    curr_hit->ddata( "time" )  : proptime;
  if ( curr_hit->mother_particle().name() == "mu+" ||
       curr_hit->mother_particle().name() == "mu-" ){
    if ( curr_hit->mother_particle().fetch_sproperty("CreatorProcess")=="none" )
      muProp++;
  } else if ( curr_hit->mother_particle().name() == "lepton_shower" )
    muProp += 0.5;
  
  // std::cout<<x<<"\t"<y<<"\t"<<z<<"\t"<<std::endl;
  
  cout<<curr_hit->x()[0]<<endl;
  cout<<curr_hit->x()[1]<<endl;
  cout<<curr_hit->x()[2]<<endl;


  //Do attenuations.
  double xE1, xE2, yE1, yE2;
  //double Xphot, Yphot;
  //Assume equal amounts of energy from both views and
  // equal energy flow in both directions along strip.
  xE1 = xE2 = yE1 = yE2 = totEng/4;
  //   proptime /= T.size();
  //   std::cout<<proptime<<std::endl;

  //double slope = OctGeom == 1 ? (_detectorY - _detectorX*tan(atan(1)/2.))/
  //(_detectorY*tan(atan(1)/2.) - _detectorX) : -1.;
  double dt, dtx, dty;
  //need to take into account drift distance to closest and furthest edge.

  xE1 = xE1 * exp(-(xedge - x)/_attLength);
  xE2 = xE2 * exp(-(3*xedge-x)/_attLength);
  yE1 = yE1 * exp(-(yedge - y)/_attLength);
  yE2 = yE2 * exp(-(3*yedge-y)/_attLength);
  
  dtx = (xedge - x) < (3*xedge - x) ?
    (xedge - x)/vlight : (3*xedge - x)/vlight;
  dty = (yedge - y) < (3*yedge - y) ?
    (yedge - y)/vlight : (3*yedge - y)/vlight;

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

#include <gdml_hit_constructor.h>

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
  cout<<"_activeLength"<<_activeLength<<endl;
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

}

gdml_hit_constructor::~gdml_hit_constructor()
{

}

void gdml_hit_constructor::reset()
{
  //Clear out map.correct??
  _voxels.clear();

}

void gdml_hit_constructor::execute(const std::vector<bhep::hit*>& hits,
				   std::vector<bhep::hit*>& rec_hit, std::vector<TH1F*>& histo_vec)
{

  rawHitsTH1F = histo_vec[0];
  clusteredHitsTH1F = histo_vec[1];
  digitizedHitsTH1F = histo_vec[2];

  std::vector<bhep::hit*> un_clustered_rec_hit;

  //First clear out map.
  reset();
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;

  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );
  
  // Calculate the z-position as the average of the z of the 4 layers of scintilating bars.
  //calculate_layerZ(sortedHits);

  //sort into voxels map.
  // parse_to_map( sortedHits );
  clustering(sortedHits);

  //Make rec_hits from vox.
  construct_hits( rec_hit );
  //construct_hits( un_clustered_rec_hit );

  // Implement the clustering algorithm to cluster hits in the plane, 
  // using which we can determine a better x,y pos then before.
  // clustering(sortedHits);

  //rec_hit = un_clustered_rec_hit;
}

void gdml_hit_constructor::clustering(const std::vector<bhep::hit*>& zSortedHits)
{
  /*
    Cluster the real hits (bar positions from hits) to produce hit positions.
    Also utilize the bar overlap to be able to give an even better position.

  */

  //Take in the filled _voxels (Sorted hits by z planes).
  // Check if two hits next to eachother (in x) means that hit activated both. (passed between)
  // Check if two hits next to eachother (in y) means that hit activated both. (passed between)
  // Return a better voxelisation?

  // std::vector<bhep::hit*> clustered_hits;  
  std::vector<bhep::hit*>::const_iterator hitIt;
  std::vector<bhep::hit*> moduleHits;
  std::vector<std::vector<bhep::hit*> > moduleHitsVector;
  std::vector<std::vector<double> > clustered_hits;
 

  // Fill vectors with hits in the same module (xyxy),sorted by barPosZ.
  for (hitIt = zSortedHits.begin();hitIt != zSortedHits.end();hitIt++)
    {
      double currZ = (*hitIt)->ddata( "barPosZ" );
      double nextZ;
      rawHitsTH1F->Fill((*hitIt)->x()[2]);
      
      if(hitIt + 1 != zSortedHits.end()){ nextZ = (*(hitIt + 1))->ddata( "barPosZ" );}
      else {nextZ = currZ + 2 * _activeLength;}

      if(fabs(currZ-nextZ) < 2 * _activeLength)
	{
	  moduleHits.push_back((*hitIt));
	}
      else//Next is to far away
	{
	  moduleHits.push_back((*hitIt));
	  moduleHitsVector.push_back(moduleHits);
	  moduleHits.clear();
	}    
    }
  // Do the actually clustering
  for(int counter = 0; counter < moduleHitsVector.size(); counter++)
    {
      if(moduleHitsVector[counter].size() != 0){clusteringXY(moduleHitsVector[counter], counter);}
    }
}

void gdml_hit_constructor::clusteringXY(const std::vector<bhep::hit*> hits, int key)
{
  std::vector<bhep::hit*> X, Y;
  double z = 0;
  
  
  for(int inCounter = 0; inCounter < hits.size(); inCounter++)
    {
      cout<<hits[inCounter]->x()[0]<<"\t"
	  <<hits[inCounter]->x()[1]<<"\t"
	  <<hits[inCounter]->x()[2]<<"\t"
	  <<hits[inCounter]->idata( "IsYBar" )<<"\t"
	  <<hits[inCounter]->idata( "barNumber" )<<"\t"
	  <<hits[inCounter]->ddata( "barPosZ" )<<"\t"
	  <<hits[inCounter]->ddata( "barPosT" )<<"\t"
	  <<endl;

      z+=hits[inCounter]->ddata( "barPosZ" );
      
      if( hits[inCounter]->idata( "IsYBar" ) == 0){X.push_back(hits[inCounter]);}
      else {Y.push_back(hits[inCounter]);}
      
    }
  
  z= z/hits.size();
  int vox_x = -1;
  int vox_y = -1;
  
  if(X.size() != 0)
    {
      vox_x = calculate_new_vox_no(X);
 
    }
  if(Y.size() != 0)
    {
      vox_y = calculate_new_vox_no(Y);
 
    }

  cout<<"vox_x "<<vox_x<<endl;
  cout<<"vox_y "<<vox_y<<endl;
  _nVoxX = 47;
  
  int vox_num = vox_x + vox_y*_nVoxX;
  cout<<"vox_num "<<vox_num<<endl;
  cout<<"z "<<z<<endl;


  //for the whole vector.

  for(int cnt = 0; cnt<hits.size();cnt++)
    { 
      if ( vox_num >= 0){_voxels[z].insert( pair<int,bhep::hit*>(vox_num,hits[cnt]) );}
    }
     clusteredHitsTH1F->Fill(z);
}

int gdml_hit_constructor::calculate_new_vox_no(std::vector<bhep::hit*> hits)
{
  /*
    Take all the hits in a module (4z planes) and calculate the voxel number from it.
    Takes in only X or only Y plane hits. Does not yet handle multiple hits per plane.

  */

  int vox_num = -1;
  int currBarNum = hits[0]->idata( "barNumber" );

  if(hits.size() == 1)
    {
      //int currBarNum = hits[0]->idata( "barNumber" );

      if(currBarNum % 2 == 0) // If barNum even then back bar
	{
	  vox_num = 2*currBarNum;
	}
      else{vox_num = 2*(currBarNum -1)+2 ;}
    }
  else
    {
      int frontBarNum = -1;
      int backBarNum = -1;

      for(int cnt = 0; cnt<hits.size(); cnt++)
	{
	  currBarNum = hits[cnt]->idata( "barNumber" );
	  if(currBarNum % 2 == 0){backBarNum = currBarNum/2;}
	  else{frontBarNum = (int) (currBarNum -1)/2;}
	}
      vox_num = 2*backBarNum + 2*frontBarNum + 1;
    }

  // cout<<"vox_num: "<<vox_num<<endl;

 return vox_num;
}

void gdml_hit_constructor::calculate_layerZ(const std::vector<bhep::hit*>& hits)
{
  // Get in z-pos sorted hits.

  //  cout<<"Find plane, transbarpos: "<<curHit.ddata( "barPosT" )<<endl;

  //cout<<"Find plane, longbarpos: "<<curHit.ddata( "barPosZ" )<<endl;

  // Find the closest bars and take an average.

  std::vector<bhep::hit*>::const_iterator hitIt;
  double previousBarPosZ = _detectorLength;
  int counter;
  double sumPos;

  for (hitIt = hits.begin();hitIt != hits.end();hitIt++){ 
    double currLongBarPosZ = (*hitIt)->ddata( "barPosZ" );

    //cout<<"In digi calculate_layerZ, currZ is: "<<currLongBarPosZ<<endl;
    //cout<<"In digi calculate_layerZ, prevZ is: "<<previousBarPosZ<<endl;

    if(previousBarPosZ == _detectorLength)
      {
	previousBarPosZ = currLongBarPosZ;
	sumPos =  currLongBarPosZ;
	counter = 1;
      }
    else
      {
	if( fabs(currLongBarPosZ - previousBarPosZ) < 4 * _activeLength)
	  {
	    sumPos +=  currLongBarPosZ;
	    counter++;
	  }
	else
	  {
	    _zLayer.push_back( sumPos/counter );
	    //cout<<"In digi calculate_layerZ, z is: "<<sumPos/counter<<endl;
	    
	    previousBarPosZ = currLongBarPosZ;
	    sumPos = currLongBarPosZ;
	    counter = 1;
	  }
      }
    
  }
  _zLayer.push_back( sumPos/counter );
  //cout<<"In digi calculate_layerZ, z is: "<<sumPos/counter<<endl;


}

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

    //cout<<"zpos in parse_to_map: "<<zpos<<endl;
    
    //get Vox number;
    voxNo = calculate_vox_no( *(*hitIt) );

    //cout<<"voxNo in parse_to_map: "<<voxNo<<endl;

    //Add voxel to map (or hit to existing voxel).
    if ( voxNo >= 0)
      _voxels[zpos].insert( pair<int,bhep::hit*>(voxNo,(*hitIt)) );
    
  }
  
}

double gdml_hit_constructor::find_plane(bhep::hit& curHit)
{
  //find the appropriate z position by comparison to
  //Layer z.
  
  double modDiff = fabs( curHit.x()[2] - (*_zIt) );

  //std::cout<<"X = "<<curHit.x()[0]<<", Y = "<<curHit.x()[1]<<", Z of hit "<<curHit.x()[2]<<std::endl;
  while ( (int)modDiff > (int) 4*_activeLength/2. && _zIt != _zLayer.end()){

    _zIt++;

    modDiff = fabs( curHit.x()[2] - (*_zIt) );
    
  }
  //cout<<"In digi find_plane, z is: "<<(*_zIt)<<endl;
  return (*_zIt);
}

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
	
	//std::cout<<"In construct_hits in hit_constructor "<<"x = "<<vhit->x()[0]<<",y = "<<vhit->x()[1]<<", z = "<<vhit->x()[2]<<std::endl;
      }
      vIt->second.erase( vIt2->first );
    }
  }
  
}

bhep::hit* gdml_hit_constructor::get_vhit(int vox, double z,
				     const std::multimap<int,bhep::hit*>& map1)
{
  //Makes a rec_hit from the voxel position and adds the relevant points.
  bhep::hit* returnPointer;  

  double totEng = 0., muProp = 0.; //these will be done in rec_hit eventually.
  vdouble X, Y, Z, E, T; //annoying but again all in rec_hit class.
  double proptime=9999999.9, vlight = 299792458. / 1.6;
  double meanvoxtime;
  int irow = vox / _nVoxX;
  int icol = vox % _nVoxX;
  
  double voxX = icol*_voxXdim + _voxXdim/2 - _detectorX/2;
  double voxY = _detectorY/2 - (irow*_voxYdim + _voxYdim/2);
  double smearingFactor = 0.06;
  double dt, dtx, dty;
  double xedge = _detectorX/2.;
  double yedge = _detectorY/2.;
  double xE1, xE2, yE1, yE2;
  double Xphot, Yphot;

  Point3D hitPos( voxX, voxY, z );

  digitizedHitsTH1F->Fill(z);
  
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
      //	       <<(*hIt).second->ddata( "EnergyDep" )<<std::endl;
    }
  //Do attenuations.

  //Assume equal amounts of energy from both views and
  // equal energy flow in both directions along strip.
  //cout<<"totEng: "<<totEng<<endl;
  xE1 = xE2 = yE1 = yE2 = totEng/2;

  //cout<<"attLength: "<<_attLength<<endl;

  xE1 = xE1 * exp(-(xedge - fabs(voxX))/_attLength);
  xE2 = xE2 * exp(-(3*xedge-fabs(voxX))/_attLength);
  yE1 = yE1 * exp(-(yedge - fabs(voxY))/_attLength);
  yE2 = yE2 * exp(-(3*yedge-fabs(voxY))/_attLength);
  
  dtx = (xedge - fabs(voxX)) < (3*xedge - fabs(voxX)) ?
    (xedge - fabs(voxX))/vlight : (3*xedge - fabs(voxX))/vlight;
  dty = (yedge - fabs(voxY)) < (3*yedge - fabs(voxY)) ?
    (yedge - fabs(voxY))/vlight : (3*yedge - fabs(voxY))/vlight;

  dt = dtx < dty ? dtx : dty;

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
      //cout<<"whit properties added"<<endl;
      
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
				     

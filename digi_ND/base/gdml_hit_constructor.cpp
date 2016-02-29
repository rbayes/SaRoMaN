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
  //_vertexDetX = store.fetch_dstore("vertex_x") * m;
  //_vertexDetY = store.fetch_dstore("vertex_y") * m;
  //_passiveLength = store.fetch_dstore("passive_thickness") * cm;
  _activeLength = store.fetch_dstore("active_thickness") * cm;
  cout<<"_activeLength"<<_activeLength<<endl;
  //_braceLength = 0.0;
  //if (store.find_dstore("bracing_thickness"))
  //  _braceLength = store.fetch_dstore("bracing_thickness") * cm;
  //_gapLength = store.fetch_dstore("air_gap") * cm;
  //_nActive = store.fetch_istore("active_layers");
  //OctGeom = 0;
  //if (store.find_istore("isOctagonal"))
  //  OctGeom = store.fetch_istore("isOctagonal");
  //_voxXdim = store.fetch_dstore("rec_boxX") * cm;
  //_voxYdim = store.fetch_dstore("rec_boxY") * cm;
  //_nVoxX = (int)( _detectorX / _voxXdim );
  //_nVox = _nVoxX * (int)( _detectorY / _voxYdim );
  _nVoxX =store.fetch_istore("nVoxX");

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
  /*
    Main executible function for the hit_constructor.
  */

  //Set the histogram plotters
  rawHitsTH1F = histo_vec[0];
  clusteredHitsTH1F = histo_vec[1];
  digitizedHitsTH1F = histo_vec[2];
  xeTH1F = histo_vec[3];
  xeAttTH1F = histo_vec[4];
  xeSmearTH1F = histo_vec[5];
  yeTH1F = histo_vec[6];
  yeAttTH1F = histo_vec[7];
  yeSmearTH1F = histo_vec[8];

  //First clear out map.
  reset();
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;
  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );

  clustering(sortedHits);

  //Make rec_hits from vox.
  construct_hits( rec_hit );
}

void gdml_hit_constructor::clustering(const std::vector<bhep::hit*>& zSortedHits)
{
  /*
    Cluster the real hits (bar positions from hits) to produce hit positions.
    Also utilize the bar overlap to be able to give an even better position.
    Main jobs is done by calling clusteringXY.
  */  
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
      else {nextZ = currZ + 3./2. * _activeLength;}

      if(fabs(currZ-nextZ) < 3./2. * _activeLength)
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
  /*
    Cluster the real hits (bar positions from hits) to produce hit positions.
    Also utilize the bar overlap to be able to give an even better position.
    Results in _voxel being filled.
  */  
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
      vox_x = calculate_vox_no(X);
 
    }
  if(Y.size() != 0)
    {
      vox_y = calculate_vox_no(Y);
 
    }

  cout<<"vox_x "<<vox_x<<endl;
  cout<<"vox_y "<<vox_y<<endl;

  int vox_num = vox_x + vox_y*_nVoxX;
  cout<<"vox_num "<<vox_num<<endl;
  cout<<"z "<<z<<endl;


  //for the whole vector.

  for(int cnt = 0; cnt<hits.size();cnt++)
    { 
      if ( vox_num >= 0){
	_voxels[z].insert( pair<int,bhep::hit*>(vox_num,hits[cnt]) );
	clusteredHitsTH1F->Fill(hits[cnt]->x()[2]);
      }
    }
    
}

int gdml_hit_constructor::calculate_vox_no(std::vector<bhep::hit*> hits)
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

  double smearingFactor = 1;
  double dt, dtx, dty;
  double xedge = _detectorX/2.;
  double yedge = _detectorY/2.;
  double xE1, xE2, yE1, yE2;
  double Xphot, Yphot;

  vector<double> barPosX;
  vector<double> barPosY;
  double sumBarPosX = 0;
  double sumBarPosY = 0;
  double barX;
  double barY;

  bhep::hit* vhit = new bhep::hit( "tracking" );

  vhit->add_property( "voxel", vox );


  // Pushback all the relevant informatiomn to the voxels. 
  // Also calculate the average position of the hits in each voxel.
  std::multimap<int,bhep::hit*>::const_iterator hIt;
  for (hIt = map1.equal_range(vox).first;hIt != map1.equal_range(vox).second;hIt++)
    {
      if((*hIt).second->idata( "IsYBar" ) == 0)
	{
	  barPosX.push_back((*hIt).second->ddata( "barPosT" ));
	  sumBarPosX+=(*hIt).second->ddata( "barPosT" );
	}
      else
	{
	  barPosY.push_back((*hIt).second->ddata( "barPosT" ));
	  sumBarPosY+=(*hIt).second->ddata( "barPosT" );
	}

      X.push_back( (*hIt).second->x()[0] );
      Y.push_back( (*hIt).second->x()[1] );
      Z.push_back( (*hIt).second->x()[2] );
      digitizedHitsTH1F->Fill((*hIt).second->x()[2]);
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
  barX=sumBarPosX/barPosX.size();
  barY=sumBarPosY/barPosY.size();

  Point3D hitPos(barX , barY, z );
  vhit->set_point( hitPos );

  //Do attenuations.

  //Assume equal amounts of energy from both views and
  // equal energy flow in both directions along strip.
  //cout<<"totEng: "<<totEng<<endl;
  xE1 = xE2 = yE1 = yE2 = totEng/2;

  xeTH1F->Fill(xE1+xE2);
  yeTH1F->Fill(yE1+yE2);

  //cout<<"attLength: "<<_attLength<<endl;

  xE1 = xE1 * exp(-(xedge + fabs(barX))/_attLength);
  xE2 = xE2 * exp(-(3*xedge-fabs(barX))/_attLength);
  yE1 = yE1 * exp(-(yedge + fabs(barY))/_attLength);
  yE2 = yE2 * exp(-(3*yedge-fabs(barY))/_attLength);

  xeAttTH1F->Fill(xE1+xE2);
  yeAttTH1F->Fill(yE1+yE2);
  
  dtx = (xedge + fabs(barX)) < (3*xedge - fabs(barX)) ?
    (xedge + fabs(barX))/vlight : (3*xedge - fabs(barX))/vlight;
  dty = (yedge + fabs(barY)) < (3*yedge - fabs(barY)) ?
    (yedge + fabs(barY))/vlight : (3*yedge - fabs(barY))/vlight;

  dt = dtx < dty ? dtx : dty;

  //double xE = xE1 + _ranGen.Gaus( 0, smearingFactor * xE1 )
  //  + xE2 + _ranGen.Gaus( 0, smearingFactor * xE2 );
  double xE = xE1 + _ranGen.PoissonD(smearingFactor * xE1)
    + xE2 + _ranGen.PoissonD(smearingFactor * xE2);
  //double yE = yE1 + _ranGen.Gaus( 0, smearingFactor * yE1 )
  //  + yE2 + _ranGen.Gaus( 0, smearingFactor * yE2 );
  double yE = yE1 + _ranGen.PoissonD(smearingFactor * yE1)
    + yE2 + _ranGen.PoissonD(smearingFactor * yE2);

  xeSmearTH1F->Fill(xE);
  yeSmearTH1F->Fill(yE);

  //if ( fabs(z) > (_detectorLength + _vertexDetdepth)/2. && 
  //     xE < _minEng && yE < _minEng )
  if ( fabs(z) > (_detectorLength)/2. && xE < _minEng && yE < _minEng )
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
				     

#include <hit_clusterer.h>

#include <CLHEP/Random/RandGauss.h>
#include <recpack/RecPackManager.h>
#include <map>

hit_clusterer::hit_clusterer(const bhep::gstore& store)
{
  //get relevant information from the store.
  // long seed = (long)store.fetch_dstore("Gen_seed");
//   _ranGen = RanluxEngine( seed, 4 );
  //_sigMa = store.fetch_dstore("pos_sig") * cm;
  _sigMaX = store.fetch_dstore("xpos_sig") * cm;
  _sigMaY = store.fetch_dstore("ypos_sig") * cm;
  _sigMaZ = store.fetch_dstore("zpos_sig") * cm;

  _res[0] = 0.85; _res[1] = 0.80; _res[2] = 0.75;

  _minEng = store.fetch_dstore("min_eng") * MeV;

  _vInX = (int)( (store.fetch_dstore("MIND_x")*m) / (store.fetch_dstore("rec_boxX")*cm) );

  _cov = EMatrix(3,3,0);
  //for (int i = 0;i < 2;i++)
  //_cov[i][i] = pow( _sigMa, 2 );
  //for (int i = 0;i < 2;i++) 
  //_cov[i][i] = pow( _sigMa, 2 );

  _cov[0][0] = pow( _sigMaX, 2 );
  _cov[1][1] = pow( _sigMaY, 2 );
  _cov[2][2] = _sigMaZ * _sigMaZ;
  
  
  _measType = store.fetch_sstore("meas_type");
  
}

hit_clusterer::~hit_clusterer()
{
}

void hit_clusterer::execute(const std::vector<bhep::hit*>& deps,
			    std::vector<cluster*>& clusts)
{
  //take a vector of voxel hits and make clusters.
  double zplane;
  std::vector<bhep::hit*>::const_iterator depIt;
  std::map<int,bhep::hit*> plane_hits;

  zplane = (*deps.begin())->x()[2];
  // std::cout<<"z of plane is "<<zplane<<std::endl;
  for (depIt = deps.begin();depIt < deps.end();depIt++){
    //std::cout<<"x = "<<(*depIt)->x()[0]<<", y = "
    //     <<(*depIt)->x()[1]<<", z = "<<(*depIt)->x()[2]<<std::endl;
    if ( (*depIt)->x()[2] != zplane ){

      //have all hits in this plane so add the clusters to vector.
      clusters_in_plane( zplane, plane_hits, clusts );

      plane_hits.clear();
      
      zplane = (*depIt)->x()[2];
      plane_hits.insert( pair<int,bhep::hit*>((*depIt)->idata("voxel"),(*depIt)) );
      
      // std::cout<<"z of plane is "<<zplane<<std::endl;
    } else plane_hits.insert( pair<int,bhep::hit*>((*depIt)->idata("voxel"),(*depIt)) );

  }
  //And the last plane.
  
  if ( plane_hits.size() != 0 ){
    clusters_in_plane( zplane, plane_hits, clusts );
    plane_hits.clear();
  }
  
}

void hit_clusterer::clusters_in_plane(double zpos, std::map<int,bhep::hit*>& deps,
				      std::vector<cluster*>& clusts)
{
  //if only a single voxel is hit in the plane make the 'cluster' directly.
  if ( deps.size() == 1 ){

    cluster* clst = make_cluster( zpos, deps.begin()->second );

    clusts.push_back( clst );

  } else {
    //Need to find out how many clusters there are and there positions.
    form_clusters( zpos, deps, clusts );
  }

}

void hit_clusterer::form_clusters(double zpos, std::map<int,bhep::hit*>& deps,
				  std::vector<cluster*>& clusts)
{
  //Have to group voxels.
  int vox1 = deps.begin()->first;
  int max_vox;

  std::map<int,bhep::hit*>::iterator mIt;

  std::map<int,double> clhit;
  clhit.insert( pair<int,double>(vox1, deps.begin()->second->fetch_dproperty("TotalEng")) );
  std::map<int,double>::iterator dIt;

  int isearch[] = {_vInX+1,_vInX,_vInX-1,1,-1,-(_vInX-1),-_vInX,-(_vInX+1)};
  
  while ( deps.size() != 0 ) {
 
    for (int ineigh = 0;ineigh < 8;ineigh++){

      mIt = deps.find( vox1 - isearch[ineigh] );

      if ( mIt != deps.end() )
	clhit.insert( pair<int,double>(vox1 - isearch[ineigh],
				       (*mIt).second->fetch_dproperty("TotalEng")) );

    }
    if ( clhit.size() == 1){
      cluster* cl1 = make_cluster( zpos, deps[clhit.begin()->first] );
      cl1->set_VoxX( 1 );
      cl1->set_VoxY( 1 );
      clusts.push_back( cl1 );
      deps.erase( clhit.begin()->first );
      clhit.clear();

      if ( deps.size() != 0 ){
	vox1 = deps.begin()->first;
	clhit.insert( pair<int,double>(vox1,
				       deps.begin()->second->fetch_dproperty("TotalEng")) );
      }

    } else {
      max_vox = get_max_vox( clhit );
      if ( max_vox == vox1 ){
	EVector pos(3,0);
	pos[2] = zpos;
	std::vector<bhep::hit*> hits;
	for (dIt = clhit.begin();dIt != clhit.end();dIt++)
	  hits.push_back( deps[dIt->first] );
	calculate_clust_pos( hits, pos );
	cluster* cl2 = make_cluster( pos, hits );
	cl2->set_VoxX(_nVoxV[0]);
	cl2->set_VoxY(_nVoxV[1]);
	clusts.push_back( cl2 );
	for (dIt = clhit.begin();dIt != clhit.end();dIt++)
	  deps.erase( dIt->first );
	clhit.clear();

	if ( deps.size() != 0 ){
	  vox1 = deps.begin()->first;
	  clhit.insert( pair<int,double>(vox1,
					 deps.begin()->second->fetch_dproperty("TotalEng")) );
	}
      } else {
	vox1 = max_vox;
	clhit.clear();
	clhit.insert( pair<int,double>( max_vox, 
					deps[max_vox]->fetch_dproperty("TotalEng")) );
      }
    }
      
  }

}

int hit_clusterer::get_max_vox(const std::map<int,double>& voxes)
{
  int max;
  double max_val;
  std::map<int,double>::const_iterator vIt2;

  max = voxes.begin()->first;
  max_val = voxes.begin()->second;

  for (vIt2 = voxes.begin();vIt2 != voxes.end();vIt2++){

    if ( vIt2->second > max_val ){
      max = vIt2->first;
      max_val = vIt2->second;
    }

  }
  
  return max;
}

cluster* hit_clusterer::make_cluster(double zpos, bhep::hit* dep)
{

  EVector hit_pos(2,0);
  hit_pos[0] = dep->x()[0];
  hit_pos[1] = dep->x()[1];
  
  EVector meas_pos(3,0);
  meas_pos[0] = hit_pos[0];
  meas_pos[1] = hit_pos[1];
  meas_pos[2] = zpos;

  cluster* me = new cluster();
  me->set_name(_measType);
  //set cov if there is low energy in one of the views.
  if ( dep->fetch_dproperty("XEng") < _minEng ) _cov[0][0] = pow( 1.0*m, 2 );
  if ( dep->fetch_dproperty("YEng") < _minEng ) _cov[1][1] = pow( 1.0*m, 2 );
  /* me->set_hv(HyperVector(hit_pos,_cov));
  me->set_name("volume", "Detector");
  me->set_position( meas_pos );*/

   //// For new recpack
  EMatrix meas_cov(3,3,0);
  me->set_hv(HyperVector(meas_pos,_cov,_measType));
  me->set_name(RP::setup_volume, "mother");
  me->set_position_hv( HyperVector(meas_pos, _cov, RP::xyz) );
  ////
  
  //for multiple track
  me->set_hv("energy", HyperVector(0,0));
  me->set_hv("MuonProp", HyperVector(0,0));
  me->set_hv("hittime", HyperVector(0,0));
  me->add_hit( dep );

  return me;
}

cluster* hit_clusterer::make_cluster(const EVector& vec,
				     const std::vector<bhep::hit*>& deps)
{

  EVector hit_pos(2,0);
  hit_pos[0] = vec[0];
  hit_pos[1] = vec[1];

  cluster* me = new cluster();
  me->set_name(_measType);
  /* me->set_hv(HyperVector(hit_pos,_cov));
  me->set_name("volume", "Detector");
  me->set_position( vec );*/

  ////
  EMatrix meas_cov(3,3,0);
  me->set_hv(HyperVector(vec,_cov,_measType));
  me->set_name(RP::setup_volume, "mother");
  me->set_position_hv( HyperVector(vec, _cov, RP::xyz) );
  ///
  
  //for multiple track
  me->set_hv("energy", HyperVector(0,0));
  me->set_hv("MuonProp", HyperVector(0,0));
  me->set_hv("hittime", HyperVector(0,0));
  std::vector<bhep::hit*>::const_iterator depIt2;
  for (depIt2 = deps.begin();depIt2 != deps.end();depIt2++)
    me->add_hit( (*depIt2) );

  return me;
}

void hit_clusterer::calculate_clust_pos(const std::vector<bhep::hit*>& hits, EVector& vec)
{
  //Weighted mean for x and y position.
  double X = 0, Y = 0, Qm1 = 0, Qm2 = 0, Q1, Q2;
  //To calculate number voxes in x and y.
  int nX = 0, nY = 0;
  std::map<double,bool> xmast, ymast;

  std::vector<bhep::hit*>::const_iterator hitIt;
  for (hitIt = hits.begin();hitIt != hits.end();hitIt++){
    
    Q1 = (*hitIt)->ddata("XEng");
    Q2 = (*hitIt)->ddata("YEng");

    if ( Q1 >= _minEng ){
      X += Q1 * (*hitIt)->x()[0];
      Qm1 += Q1;
    }
    if ( Q2 >= _minEng ){
      Y += Q2 * (*hitIt)->x()[1];
      Qm2 += Q2;
    }

    if ( hitIt == hits.begin() ){
      xmast[(*hitIt)->x()[0]] = true;
      if ( (*hitIt)->fetch_dproperty("XEng") >= _minEng ) nX++;
      ymast[(*hitIt)->x()[1]] = true;
      if ( (*hitIt)->fetch_dproperty("YEng") >= _minEng ) nY++;
    } else {
      if ( xmast.find( (*hitIt)->x()[0] ) != xmast.end() &&
	   (*hitIt)->fetch_dproperty("XEng") >= _minEng ){
	xmast[(*hitIt)->x()[0]] = true; nX++;}

      if ( ymast.find( (*hitIt)->x()[1] ) != ymast.end() &&
	   (*hitIt)->fetch_dproperty("YEng") >= _minEng ){
	ymast[(*hitIt)->x()[1]] = true; nY++;}
    }

  }
  _nVoxV[0] = nX; _nVoxV[1] = nY;
  if ( nX == 0 ) _cov[0][0] = pow( 1.0*m, 2 );
  else if ( nX > 1 ) _cov[0][0] = pow( _res[nX-1]*_sigMaX, 2 );
  if ( nY == 0 ) _cov[1][1] = pow( 1.0*m, 2 );
  else if ( nY > 1 ) _cov[1][1] = pow( _res[nY-1]*_sigMaY, 2 );

  vec[0] = X / Qm1;// + RandGauss::shoot(&_ranGen, 0, _sigMa);
  vec[1] = Y / Qm2;// + RandGauss::shoot(&_ranGen, 0, _sigMa);
  
  xmast.clear(); ymast.clear();
}

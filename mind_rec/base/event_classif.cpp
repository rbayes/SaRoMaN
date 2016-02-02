
/***********************************************************************
 * Implementation of neutrino event classification class.              *
 * The functions here will take in the measurements associated to an   *
 * event and attempt to reject neutral current and electron nu events  *
 * while identify CC interaction types and performing appropriate muon *
 * extraction.                                                         *
 *                                                                     *
 * Execute will return an integer code based on what type the event    *
 * has been classified as: 0 - NC/e, 1 - CC DIS, 2 - CC Quasi etc.     *
 *                                                                     *
 * Author: Andrew Laing, Feb. 2009                                     *
 ***********************************************************************/

#include <event_classif.h>
#include <recpack/stc_tools.h>
#include <TMath.h>

using namespace bhep;

class sortTrajByHits{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->size() > t2->size()) return true;
    return false;
  }
};

class sortTrajByLength{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->length() > t2->length()) return true;
    return false;
  }
};

class sortHitsByR{
public:
  bool  operator()(cluster* hit1, cluster* hit2 ){
    double r1= sqrt(pow( hit1->position()[0],2) + pow( hit1->position()[1],2));
    double r2= sqrt(pow( hit2->position()[0],2) + pow( hit2->position()[1],2));
    if (r1 > r2) return true;
    return false;
  }
};

class sortHitsByZ{
public:
  bool  operator()(cluster* hit1, cluster* hit2 ){
    if (hit2->position()[2] > hit1->position()[2]) return true;
    return false;
  } 
};

//**********************************************************************
event_classif::event_classif() {
  //**********************************************************************
}

//***********************************************************************
event_classif::~event_classif() {
  //***********************************************************************
}

//Initialization of the class
//***********************************************************************
void event_classif::Initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, double wFe, MINDsetup* geom) {
  //***********************************************************************

  _m = bhep::messenger( vlevel );
  _m.message("++++ Classifier  init  function ++++",bhep::VERBOSE);

  _infoStore = pstore;
  
  /// read parameters from file
  readParam();
  
  _voxEdge = _infoStore.fetch_dstore("rec_boxX") * cm;

  _FeWeight = wFe;

  _geom = *geom;

  cout<<"max_z: "<<_geom.getZMax()<<endl;
  cout<<"min_z: "<<_geom.getZMin()<<endl;
  
  double vertdepth = _infoStore.find_dstore("vertexDepth")?
    _infoStore.fetch_dstore("vertexDepth") : 0.0;
  _detLength = (_infoStore.fetch_dstore("MIND_z") + vertdepth)/2. * m;
  
  if( _infoStore.find_istore("max_N_trajectories") )
    _maxNtraj = _infoStore.fetch_istore("max_N_trajectories");
  else
    _maxNtraj = 8;

  if ( _outLike ){
    string likeFile = pstore.fetch_sstore("like_file");
    _outFileEv = new TFile(likeFile.c_str(), "recreate");
    _likeTree = new TTree("h1", "Possible Likelihood parameters");
    
    /// branches for likelihhod file 
    set_branches(); 
  }
}

//Execution of the class
//***********************************************************************
bool event_classif::Execute(const vector<cluster*>& hits,
			    vector<Trajectory*> &vtrajs, vector<cluster*>& hads) {
  //***********************************************************************
  
  _m.message("++++ Classifier Execute Function ++++", bhep::VERBOSE);
  
  ///reset variable for each event
  reset();
  
  ///copy hits
  std::vector<cluster*> hits2 = hits;
  cout<<"Num hits: "<<hits2.size()<<endl;

  ///start looking for trajectories
  bool ok;
  //while (hits2.size() > 5 && (int)vtrajs.size() <= _maxNtraj) {  
  while (hits2.size() >= 4 && (int)vtrajs.size() <= _maxNtraj) {  
    
    ///create trajectory     
    Trajectory* traj = new Trajectory();
    
    //reset each trajectory informations    
    reset_traj_info();
    
    ///sort hits in increasing Z positions
    sort( hits2.begin(), hits2.end(), sortHitsByZ() );

    //muontraj.set_quality("lowPt",0);
    //traj->set_quality("lowPt",0);


    cout<<"hits2: "<< hits2[0]->position()[2]<<endl;;
    cout<<"hits2: "<< hits2[hits2.size()-1]->position()[2]<<endl;
    cout<<"hits2: "<<hits2.size()<<endl;

    ///calculate the number of planes containing single hit and arrange the z-position, energy of the hits in increasing z 
    //Occupancy.
    ok = get_plane_occupancy( hits2 );

    ///called if liklihood output has to generate
    if ( _outLike )
      output_liklihood_info( hits2 );

    /// CA and PR both
    ok = chargeCurrent_analysis(hits2, *traj, hads);

    cout<<"ChargeCurrent ok: "<<ok<<endl;
    // ok = muon_extraction_through_PatterRec(hits2, *traj, hads);

    _m.message("nmeas in traj=",traj->size(),"  and hadrons left=",hads.size()," failType=",_failType, bhep::VERBOSE);
    
    ///if CA and free section KF filtering not fails   
    //if(_failType != 4 && _failType != 5 && traj->nmeas()>=5 )  {
    if(ok) 
      {
	///for liklhood
	if (_outLike)
	  traj_like( hits2 );
	
	if ( _outLike ) out_like();
	
	//fill the vector of trajectories
	vtrajs.push_back(traj);
	///set information of trajectory
	fill_traj_info( *traj); 
	// reset each trajectory informations    
	reset_traj_info();
	///clear vector of clusters
	hits2.clear();
	///if hadrons > 5 then loop
	if((int)hads.size() >= _min_hits || (int)vtrajs.size() > _maxNtraj) hits2 = hads;
	else break;
	
	///clear hadron container
	hads.clear();
      }
    else {
      hads=hits2;
      _m.message("hadron left =",hads.size(), bhep::VERBOSE);
      
      hits2.clear();
      delete traj;
      break;
    }
  }
  _m.message("eventclass::PR, size = ", _vPR_seed.size()," vtraj size =",vtrajs.size(), bhep::VERBOSE);
  
  //sort the trajectories (not working because length =0 ??)
  if(vtrajs.size()>1)  sort( vtrajs.begin(), vtrajs.end(), sortTrajByLength());

  return ok;
}

//***********************************************************************
void event_classif::Finalize() {
  //***********************************************************************
  
  if ( _outLike ){
    _outFileEv->Write();
    _outFileEv->Close();

    delete _outFileEv;
  }
}

//***********************************************************************
void event_classif::readParam() {
  //***********************************************************************
 
  _m.message("++++ readParam function of classifier ++++",bhep::VERBOSE);
  
 if ( _infoStore.find_istore("likeli") )
    _outLike = _infoStore.fetch_istore("likeli");
  else _outLike = false;
  
  _max_consec_missed_planes = _infoStore.fetch_istore("max_consec_missed_planes");
  _min_plane_prop = _infoStore.fetch_dstore("min_use_prop");
  _min_seed_hits = _infoStore.fetch_istore("min_seed_hits");
  _min_check =  _infoStore.fetch_istore("min_check_nodes");
  _min_hits = _infoStore.fetch_istore("low_Pass_hits");
  _max_hits = _infoStore.fetch_istore("plane_occupancy");  
  _max_nmult = _infoStore.fetch_istore("max_multOcc_plane");
  _tolerance = _infoStore.fetch_dstore("pos_res") * cm;
  _chi2_max = _infoStore.fetch_dstore("accept_chi");
  _max_coincedence = _infoStore.fetch_dstore("max_coincidence");
  
  _pieceLength = _infoStore.fetch_dstore("passive_thickness") * cm
    + _infoStore.fetch_dstore("active_thickness") * _infoStore.fetch_istore("active_layers") * cm
    + _infoStore.fetch_dstore("air_gap") * (_infoStore.fetch_istore("active_layers")+1) * cm
    + _infoStore.fetch_dstore("bracing_thickness") * _infoStore.fetch_istore("active_layers") * cm;
  
  _isTASD = _infoStore.fetch_dstore("passive_thickness") == 0 ? true : false;

  _maxBlobSkip = _infoStore.fetch_dstore("maxBlobSkip");
  _minBlobOcc = _infoStore.fetch_dstore("minBlobOcc");

  _detX = _infoStore.fetch_dstore("MIND_x") * m;
  _detY = _infoStore.fetch_dstore("MIND_y") * m;
  _detZ = _infoStore.fetch_dstore("MIND_z") * m;
  _vdetZ = _infoStore.find_dstore("vertex_z")?
    _infoStore.fetch_dstore("vertex_z") : 0.0;

  if(_infoStore.find_dstore("WLSatten"))
    _WLSAtten = _infoStore.fetch_dstore("WLSatten");
  else						
    _WLSAtten = 5000.;

  _octGeom  = _infoStore.find_istore("IsOctagonal") ?
    _infoStore.fetch_istore("IsOctagonal"): 1;
}

//***********************************************************************
void event_classif::reset(){
  //***********************************************************************
  //Resets relevant members at the start of each execute.

  _lowPt = 0;

  _nplanes = 0;
  _freeplanes=0;
  _meanOcc = 0;
  _badplanes = 0;
  _longestSingle = 0;
  _endProj = false;
  _vRadCount.clear();
  _radialLongest = 0;
  _Xtent = 0;
  _intType = 0;
  _endLongSing = 0;
  _endLongPlane = 0;
  _lastIso = 0;
  _vertGuess = 0;
  _vertGuessZ = 0;
  _failType = 0;

  _planes.clear();
  _radPlanes.clear();

  stc_tools::destroy(_planes);
  stc_tools::destroy(_radPlanes);
  
  _vPR_seed.clear();
  _nRadPlanes=0;
  _planeEnd.reset();
}

//***********************************************************************
bool event_classif::get_plane_occupancy(vector<cluster*>& hits){
  //***********************************************************************
  //Gets plane occupancies and total plane energies.
  //Needs hits in increasing z order.
  
  //   std::cout<<"+++++I am in get_plane_occupancy"<<std::endl;
  _m.message("++++ Calculating plane energies and occupancies ++++",bhep::VERBOSE);
  
  bool ok = true;

  /// size of vector<cluster> hits
  size_t nHits = hits.size();
  //cout<<"get_plane_occ ()  :: nHits="<<nHits<<endl;
  int single_count = 0;
  double testZ, testX, testY, curZ;
  size_t hits_used = 0, imeas = 0;
  int planeIndex=0 ;

  _longestSingle = 0;

  if ( !_outLike ) _visEng = 0;

  /// loop over hits to calculate plane occupancy  
  do {
    
    testX = hits[imeas]->position()[0];
    testY = hits[imeas]->position()[1];
    testZ = hits[imeas]->position()[2];
    hits_used++;
    // Avoid a hit with an undefined z position
    if ( fabs(testZ) >  _detLength ) 
      continue;
    // If nan.
    if(testX != testX) continue;
    if(testY != testY) continue;
    if(testZ != testZ) continue;
    
    ///create plane info
    plane_info* plane = new plane_info(planeIndex, testZ, _infoStore);
    plane->AddHit(hits[imeas]);
    cout<<"testZ: "<<testZ<<endl;
    
    ///calculate the z position which is the current z for hits 1 -> total no of hits in the cluster 
    for (size_t i = hits_used;i <nHits;i++) {
      curZ = hits[i]->position()[2];

      cout<<"curZ: "<<curZ<<" testZ: "<<testZ<<" _tolerance: "<<_tolerance<<endl;
   
      if (curZ <= testZ + _tolerance) {
	
	// add the hit to the same plane
	plane->AddHit(hits[i]);
	cout<<"Added hit"<<endl;
	
	testZ = hits[i]->position()[2];
	hits_used++;
      } else break;
      
    }
    
   
    _m.message(" get plane info =",plane->GetZ()," Occ=",plane->GetNHits()," PlaneNo=",plane->GetPlaneNo(), bhep::VERBOSE);
    
    cout<<"get plane info = "<<plane->GetZ()<<" Occ= "<<plane->GetNHits()<<" PlaneNo= "<<plane->GetPlaneNo()<<endl;
    ///fill the plane_info vector
    _planes.push_back(plane);
    
   
    ///increase the plane index
    planeIndex++;
    
    //------------------------------------------------------------     
    ///look for the largest free section position
    if ( plane->GetNHits() == 1 )
      single_count++;
    else {
      if ( single_count >= _longestSingle ){
	_longestSingle = single_count;
	_endLongSing = imeas - 1;
	_endLongPlane = (int)_planes.size()-2;
      }
      single_count = 0;
    }
    
    
    if ( imeas == nHits -1 && single_count != 0 &&  plane->GetNHits() == 1 ){
      if ( single_count >= _longestSingle ){
	_longestSingle = single_count;
	_endLongSing = imeas;
	_endLongPlane = (int)_planes.size()-1;
	
      }
    }
    //-------------------------------------------------------------------------------
    
    //total energy 
    if ( !_outLike ) _visEng += plane->GetEng();
    
    imeas +=  plane->GetNHits();
    _meanOcc += (double) plane->GetNHits();

    
  } while (hits_used != nHits);
  
  
  ///total no of planes
  _nplanes = (int)_planes.size();
  

  /// Mean occupancy
  if ( _nplanes == 0 ) return false;
  _meanOcc /= (double)_nplanes;
  
  ///if no of single-hits in the longest part is > 5 and free part is 50% of the hitsperplane then 
  ///assess the event‚¶ free section to look for viable muon
  
   
  /*//----------------------------------------------
    for(int pl = _nplanes -1; pl>=0; pl--){
    _longestSingle++;
    if(_planes[pl]->GetNHits()>1) break;
    }
  
    /// to calculate end point of the free section
    int free_count=0, count=0;
    for(int pl = _nplanes -1; pl>=0; pl--){
    
    if(free_count > 0) break;
    
    if(_planes[pl]->GetNHits()==1) {
    _endLongPlane = pl;
    _endLongSing = (nHits -1) - count;
    free_count++;
    } 
    else count++;
    }
    //cout<<"endLongPlane="<<_endLongPlane<<"  endLongSing="<<_endLongSing<<endl;
    
    //----------------------------------------------*/
  
  // if ( _longestSingle >= _min_seed_hits && _endLongPlane > 0.5*(double)_nplanes )
  if ( _longestSingle >= _min_seed_hits)
    assess_event( hits );
   
   
  return ok;
}

//***********************************************************************
void event_classif::assess_event(vector<cluster*>& hits)
  //***********************************************************************
{
   _m.message("++++ Assess event for complicated end point ++++",bhep::VERBOSE);
  
  //Assess the event 'free' sections to look for a viable
  //muon and tidy up the endpoint if neccessary.
  

  double dispTrans[2], endMeanOcc = 0, skipSize = 0;
  unsigned int evEnd = hits.size() - 1;
  
  vector<int>::iterator aIt;
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";
  
  
  ///At the last hit of the longest free section  
  if ( _endLongSing == (int)evEnd ){

    cout<< hits[evEnd]->vector()[0]<<" "<< hits[evEnd]->vector()[1]<<" "<< hits[evEnd]->vector()[2]<<endl;
    cout<< hits[evEnd-1]->vector()[0]<<" "<< hits[evEnd-1]->vector()[1]<<" "<< hits[evEnd-1]->vector()[2]<<endl;
    cout<< hits[evEnd-2]->vector()[0]<<" "<< hits[evEnd-2]->vector()[1]<<" "<< hits[evEnd-2]->vector()[2]<<endl;
    
    
    dispTrans[0] = sqrt( pow( hits[evEnd]->vector()[0]-hits[evEnd-1]->vector()[0], 2)
			 + pow( hits[evEnd]->vector()[1]-hits[evEnd-1]->vector()[1], 2) );
    dispTrans[1] = sqrt( pow( hits[evEnd-1]->vector()[0]-hits[evEnd-2]->vector()[0], 2)
			 + pow( hits[evEnd-1]->vector()[1]-hits[evEnd-2]->vector()[1], 2) );
    
    
    cout<<"dispTrans[0]: "<<dispTrans[0]<<endl;
    cout<<"dispTrans[1]: "<<dispTrans[1]<<endl;
    cout<<"_voxEdge*10: "<<_voxEdge*10<<endl;

    //Just in case there is a bad hit at the endpoint.
    
    if ( dispTrans[0] > _voxEdge*10 && dispTrans[1] <= _voxEdge*10 ){
      hits[evEnd]->set_name(skipped,hit_in);
      _badplanes = 1;
      /////
      _planeEnd = *_planes[_nplanes - 2];
      
    } else if ( dispTrans[1] > _voxEdge*10 ){
      hits[evEnd]->set_name(skipped,hit_in);
      hits[evEnd-1]->set_name(skipped,hit_in);
      _badplanes = 2;
           
      /////
      _planeEnd = *_planes[_nplanes - 3];
      
    } else {
      
      _planeEnd = *_planes[_nplanes - 1];
      
    }
  } else {
    //A more complicated endpoint.
    
    ///calculate the occupancy of the end planes
    for ( int pl =_planes.size() -1; pl >= _endLongPlane; pl--){
      endMeanOcc += (double)_planes[pl]->GetNHits();
      skipSize++;
      
    }
    endMeanOcc /= (double)skipSize;
   

    /// if occpancy of the end planes are larger than 2 (minBlobOcc)  
    if ( endMeanOcc > _minBlobOcc ){
      
      
      /// if 20% of planes are skipped OR no of free planes >=10
      if ( skipSize/(double)_planes.size() < _maxBlobSkip
	   || _longestSingle >= 2*_min_seed_hits ){
	
	
	///skip the hits of the end planes	
	for ( int pl =_planes.size() -1; pl >= _endLongPlane; pl--){
	  for(int ht=0; ht<_planes[pl]->GetNHits(); ht++){
	    (_planes[pl]->GetHits()[ht])->set_name(skipped,hit_in);
	  }
	}
	///
	_planeEnd = *_planes[_endLongPlane];
	
	_badplanes = (int)skipSize;
      } else _longestSingle = 0; //Brute force way to reject.
    } else {
      
      _endProj = true;
      
      _planeEnd = *_planes[_endLongPlane];
    }
    
  }
  _m.message("planeEnd info :: end planeNo=",_planeEnd.GetPlaneNo(), bhep::VERBOSE);
}

//***********************************************************************
bool event_classif::chargeCurrent_analysis(vector<cluster*>& hits,
						      Trajectory& muontraj, vector<cluster*>& hads){
  //***********************************************************************
  _m.message("++++ Performing muon reconstruction ++++",bhep::VERBOSE);
  //cout<<"inside muon_extraction_through_PatterRec"<<endl;
  
  muontraj.reset();
  bool ok = true;
  
  _recChi = EVector(3,0);
  _recChi[1] = 100000; //SetLarge dummy value for minChiHadron bit.
  
  /*//looking into beginning of planes to guess vertex hitno, to exclude backward particle hits and 
if not found then  excluded_hits = 0; _exclPlanes = 0; i.e, vertGuess =0*/
  //_vertGuess = exclude_backwards_particle();
  
  ///Zpos of vertex
  //_vertGuessZ = hits[_vertGuess]->position()[2];

  //Zpos of the start of hits.
  _vertGuessZ = hits[0]->position()[2];

  cout<<"_vertGuessZ: "<<_vertGuessZ<<endl;
  
  //_m.message("_vertGuess =",_vertGuess,"  Z= ",hits[_vertGuess]->position()[2], bhep::VERBOSE);
  
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  const dict::Key not_in = "False";
  
  
  //if ( _longestSingle >= _min_seed_hits && _endLongPlane > 0.5*(double)_nplanes ){
  ///add single occupancy planes measurements to muontraj 
  //for(int pl=_planeEnd.GetPlaneNo(); pl >= _exclPlanes+_min_check; pl-- ){
  cout<<"_planeEnd: "<<_planeEnd.GetPlaneNo()<<" size: "<<_planes.size()<<endl;
  //for(int pl=_planeEnd.GetPlaneNo(); pl >= _min_check; pl-- ){
  for(int pl=_planes.size()-1; pl >= 0; pl-- ){
    cout<<"Planes z: "<<_planes[pl]->GetZ()<<endl;
    
    if(_planes[pl]->GetNHits()!=1) break;
    
    RecObject* ro = dynamic_cast<RecObject*>(&(*(_planes[pl]->GetHits()[0])));
    muontraj.add_node(Node(*ro));
    (_planes[pl]->GetHits()[0])->set_name(candHit, hit_in);
  }  
  
  _m.message("The free section muontraj=",muontraj, bhep::DETAILED);

  cout<<"The free section muontraj= "<<muontraj.size()<<endl; 
  
  //information from the trajectory 
  _lastIso = (int)muontraj.size();
  _m.message(" _lastIso = ",_lastIso,bhep::DETAILED);
  
  ///freeplanes inside the candidate muon trajectory
  if ( !_outLike ) _freeplanes = (int)muontraj.size();
  _m.message(" _freeplanes =",_freeplanes,"  _meanOcc=",_meanOcc, bhep::DETAILED);
  
  //set to reconstruction mode.
  if ( !MINDfitman::instance().in_rec_mode() )
    MINDfitman::instance().rec_mode();
  
  //for free tracks, _meanOcc ==1;
  if ( _meanOcc == 1)
    {
      _intType = 2;
      // Do we have enough hits to perform a curve fit?
      //if (_longestSingle >= _min_seed_hits){
      //if((int)muontraj.size() >= _min_seed_hits){ok = muon_extraction( hits, muontraj, hads);} 
      // 6Hits is the minimal needed to fit a helix.
      if((int)muontraj.size() >= 6){ok = muon_extraction( hits, muontraj, hads);} 
      //Low momentum track, we need atleast 4 hits to be able to get momenta and charge.
      else if((int)muontraj.size() >= 4){ok = LowMomentumExtraction( hits, muontraj, hads);}
      else ok = false;
      
    } 
  else 
    {
      if ( (int)muontraj.size() < _min_seed_hits ) {
	_intType = 5;
	if (!_isTASD)
	  ok = invoke_cell_auto( hits, muontraj, hads);
	else ok = false;
	if ( !ok ) _failType = 4;
	
	if(!ok)  _m.message(" invoke_cell_auto not ok", bhep::VERBOSE);
	
      } 
      else {
	if ( _badplanes !=0 ) _intType = 3;
	ok = muon_extraction( hits, muontraj, hads);
	
      }
      
    }
   
  return ok;
}

//***********************************************************************
bool event_classif::muon_extraction(vector<cluster*>& hits,
				    Trajectory& muontraj, vector<cluster*>& hads) {
  //***********************************************************************
  
  //Decide on appropriate pattern recognition and perform it.
  //Are there ways to avoid doing full extraction algorithm.
  
  _m.message("++++ event_classif::muon_extraction +++++++++++++", bhep::VERBOSE);
 
  bool ok;

  ///if vertex is otherthan 0th hit position, then those hits will be candidate hadron 
  if (_vertGuess != 0)
    for (int i = 0;i < _vertGuess;i++){
      const dict::Key candHit = "inhad";
      const dict::Key hit_in = "True";
      hits[i]->set_name(candHit, hit_in);
      hads.push_back( hits[i] );
    
      _m.message(" 1:hits[hits.size()]->position()[2] ",hits[hits.size()-1]->position()[2],"  hits[hits.size()-1]->position()[2]=",hits[hits.size()-2]->position()[2], bhep::DETAILED);
    }
  
  State patternSeed;
  
  ok = get_patternRec_seed( patternSeed, muontraj, hits);
  if(!ok)  _m.message(" get_patternRec_seed not ok",bhep::DETAILED); 
  
  if ( ok )
    ok = perform_muon_extraction( patternSeed, hits, muontraj, hads);
  if(!ok) _m.message("perform_muon_extraction not ok",bhep::DETAILED);

  ///assign the seed state
  if ( ok )
    _seedState = patternSeed;
  
  if(ok) _m.message(" event_class: traj nmeas=",muontraj.size()," intType =",_intType,"  && PR seed is=",_seedState,bhep::DETAILED);


  /// to get all the PR seed for reseed_traj inside fitter
  
  if(ok) _vPR_seed.push_back(_seedState);
  
  
  return ok;
}

//***********************************************************************
bool event_classif::get_patternRec_seed(State& seed, Trajectory& muontraj,
					vector<cluster*>& hits) {
  //***********************************************************************
  //std::cout<<"++++ event_classif::get_patternRec_seed +++++++++++++"<<std::endl;

  
  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);
  
  ///x,y,z values of the 1st node
  V[0] = muontraj.nodes()[0]->measurement().vector()[0];
  V[1] = muontraj.nodes()[0]->measurement().vector()[1];
  V[2] = muontraj.nodes()[0]->measurement().position()[2];

  //direction
  double dqtot = fit_parabola( V, muontraj);
  
  //Momentum. Estimate from empirical extent function.
  ///if cluster contains hits then _Xtent is the diff in z b/w 0th and last hit in cluster 
  ///if not, then _Xtent is the diff in z b/w 1st and one before last node (? not last node)

  /*
  if ( hits.size() != 0 ){

    _Xtent = hits[_endLongSing]->position()[2]
      - hits[_vertGuess]->position()[2];
    
  } else {
     
    _Xtent = muontraj.nodes()[(int)muontraj.size()-1]->measurement().position()[2]
      - muontraj.nodes()[0]->measurement().position()[2]; 
  }
  
  _m.message("xtent=",_Xtent, bhep::VERBOSE);
  */
  //2 pSeed = 668 + 1.06*_Xtent; //estimate in MeV, log for fit.
  // double pSeed = (9180-6610*_FeWeight) + (-2.76+4.01*_FeWeight)*_Xtent; //best for 3/2 seg


  //Use slope in bending plane to try and get correct sign for momentum seed.
  // pSeed *= -fabs(V[3])/V[3];//Pos. gradient in X is ~negative curvature.
  // Assume that the R-Z plane is the bending plane for the toroidal field
  double VRad =  dqtot/fabs(dqtot);
  //if( VRad != 0 ) 
  //pSeed *= fabs(VRad)/VRad;
  
  //V[5] = 1./pSeed;

  cout<<"In get_patternRec_seed: "<<V[5]<<" "<<1/V[5]<<endl;
  
  //Errors
  M[0][0] = M[1][1] = 15.*cm*cm;
  M[2][2] = EGeo::zero_cov()/2;
  M[3][3] = M[4][4] = 1.5;
  M[5][5] = pow(V[5],2)*4;

  //Sense
  V2[0] = 1;
  
  //Seedstate fit properties
  seed.set_name(RP::particle_helix);
  seed.set_name(RP::representation,RP::slopes_curv_z); 
  seed.set_hv(RP::sense,HyperVector(V2,M2,RP::x));
  seed.set_hv(HyperVector(V,M,RP::slopes_curv_z));
  
  // std::cout<<seed.hv().representation()<<std::endl;

  // std::cout<<"Is good here"<<std::endl;
  //  man().model_svc().model(RP::particle_helix).representation(RP::slopes_curv_z)
  //     .convert(seed,RP::default_rep);
  
  bool ok = perform_kalman_fit( seed, muontraj);
  
  if ( !ok )
    _failType = 5;
  
  return ok;
}





//***********************************************************************
double event_classif::fit_parabola(EVector& vec, Trajectory& track) {
  //***********************************************************************
  // No longer well named

  int fitcatcher;
  int nMeas = track.size();

  //if (nMeas > 4) nMeas = 4;

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];/// r[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;
  //double pathlength=0;/// tminR=9999.9;

  double firstNodeZ = track.nodes()[nMeas-1]->measurement().position()[2];
  cout<<"firstNodeZ: "<<firstNodeZ<<endl;

  double pathlength=track.nodes()[0]->measurement().position()[2] - firstNodeZ;
  
  
  pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
  pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
  z[nMeas-1] = track.nodes()[nMeas-1]->measurement().position()[2];
  pos[2] = 0.0;

  for (int iMeas = nMeas-2;iMeas >= 0;iMeas--){
    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];
    // get the b-field from the previous step
    EVector B = _geom.getBField(pos);
    u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
    pos[0] = x[iMeas];  pos[1] = y[iMeas];   pos[2] = z[iMeas];
    EVector dR = EVector(3,0);
    dR[0] = x[iMeas+1] - x[iMeas];
    dR[1] = y[iMeas+1] - y[iMeas];
    dR[2] = z[iMeas+1] - z[iMeas];
    double dq = 0.;
    //pathlength += dR.norm();
    if(iMeas < nMeas-3){
      EVector dR1 = EVector(3,0);
      dR1[0] = x[iMeas+2] - x[iMeas+1];
      dR1[1] = y[iMeas+2] - y[iMeas+1];
      dR1[2] = z[iMeas+2] - z[iMeas+1];
      EVector ddR = dR1 + dR;
      dq = dot(ddR, crossprod(Z,B))/ (crossprod(Z,B).norm());
    }
    if(pdR != 0)
      if(// minR < R && dR/fabs(dR) == -pdR/fabs(pdR) && 
	 (x[iMeas]/fabs(x[iMeas]) != x[iMeas+1]/fabs(x[iMeas+1]) ||
	  y[iMeas]/fabs(y[iMeas]) != y[iMeas+1]/fabs(y[iMeas+1]))){ 
	// Sign change and minimum momentum
	minR = pos.norm();
	minindex = iMeas;
	sumdq = 0.0;
      }
    pdR = dR.norm();
    sumdq += dq;
  }
  //cout<<"sumdq/(nMeas-2): "<<sumdq/(nMeas-2)<<endl;
  
  
  ///double wFe = geom.get_Fe_prop();
  //double p = RangeMomentum(pathlength);

  //double p2 = MomentumFromDeflection(track,0);

  double p = RangeMomentum(pathlength,firstNodeZ);

  

  //TGraph *gr1 = new TGraph((const int)minindex, z, x);
  //TGraph *gr2 = new TGraph((const int)minindex, z, y);
  //TGraph *gr3 = new TGraph((const int)minindex, z, u);

  //TF1 *fun = new TF1("parfit","[0]+[1]*x",-3,3);
  //fun->SetParameters(0.,0.001);
  TF1 *fun2 = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun2->SetParameters(0.,0.001,0.001);
  
  /*fitcatcher = gr1->Fit("parfit", "QN");
  vec[3] = fun->GetParameter(1);

  fun->SetParameters(0.,0.001);
  fitcatcher = gr2->Fit("parfit", "QN");
  vec[4] = fun->GetParameter(1);

  // fun->SetParameters(0.,0.001);
  */
  //fitcatcher = gr3->Fit("parfit2", "QN");
  double qtilde = -0.3*pow(1+fun2->GetParameter(1),3./2.)/
    (2*fun2->GetParameter(2));
  
  // int meansign = sumdq/fabs(sumdq);
  int meansign = (int)(qtilde/fabs(qtilde));

  vec[5] = meansign/p;

  //delete gr1;
  //delete gr2;
  //delete fun;
  // return sumdq;
  return qtilde;

}

//***********************************************************************
void event_classif::set_de_dx(double mom){
  //***********************************************************************
  
  double de_dx = -( 12.37 * _FeWeight * pow( mom, 0.099) 
		    + (1 - _FeWeight) * 2.16 * pow( mom, 0.075) );
  de_dx *= MeV/cm;
  
  man().geometry_svc().setup().set_volume_property_to_sons("mother","de_dx",de_dx);

}




//***********************************************************************
bool event_classif::perform_kalman_fit(State& seed, Trajectory& track) {
  //***********************************************************************
  _m.message(" event_classif :: perform_kalman_fit", bhep::VERBOSE);    
  


  ///fit the track using the seed state                             
  bool ok = man().fitting_svc().fit(seed, track);
 
  ///print first state here, which is the 1st measurement in the trajectory
  if (ok)
    seed = track.state(track.first_fitted_node());

  ///if ok then perform_muon_extraction, if not ok failtype=5
  return ok;
  
}





//***********************************************************************
bool event_classif::perform_muon_extraction(const State& seed, vector<cluster*>& hits,
					    Trajectory& muontraj, vector<cluster*>& hads) {
  //***********************************************************************


  //Loop through multiple occupancy planes finding the best match to the muon
  //traj and filtering it into the trajectory.

  _m.message("I am  event_classif::perform_muon_extraction ", bhep::DETAILED); 


  bool ok;
  long ChiMin;
  int nConsecHole = 0;
 
   
  /// start adding mulp. occupancy planes hits from the LAST PLANE of "free section" muontraj
  
  int plane_start = _nplanes -  muontraj.size() -1;
  for (int pl = plane_start; pl>=0; pl--){
    
    _m.message("plane no=",pl,"  & hits = ",_planes[pl]->GetNHits()," & intType = ",_intType, bhep::DETAILED);
    
    /// if hit corresponds to vertex position
    if (_planes[pl]->GetZ() < hits[_vertGuess]->position()[2]) break;
    
    double Chi2[(const int)(_planes[pl]->GetNHits())] ;
    
    
    ///loop over hits on a plane
    for (int ht=0; ht<_planes[pl]->GetNHits(); ht++){
      // _m.message("meas ==> ",(_planes[pl]->GetHits()[ht]->vector()),bhep::VERBOSE); 
      try {
	ok = man().matching_svc().match_trajectory_measurement(muontraj, *(_planes[pl]->GetHits()[ht]), Chi2[ht]);
      } catch (const char* msg){
	ok = false;
	std::cout<<msg<<std::endl;
      }
      // _m.message("meas.dim() = ",(_planes[pl]->GetHits()[ht])->dim(),bhep::VERBOSE);
      if ( !ok ) Chi2[ht] = 9999999;
      
    }
    
    
    /// Return hit index having min chi2
    ChiMin = TMath::LocMin((const int)(_planes[pl]->GetNHits()), Chi2);
    
    
    
    ///The hit corresponds to min Chi2 will be added to the Trajectory & assigns as muon candidate, else considered as Hadron
    for (int iht = 0; iht < _planes[pl]->GetNHits();iht++){
      
      if (iht==(int)ChiMin) {
	
	ok = man().fitting_svc().filter(*(_planes[pl]->GetHits()[iht]), seed, muontraj);
	_m.message("plane No =",_planes[pl]->GetPlaneNo()," & the hit with Z=",_planes[pl]->GetHits()[iht]->position()[2], bhep::DETAILED);
	
	///assign as muon candidate	    
	if ( ok ) {
	  const dict::Key candHit = "inMu";
	  const dict::Key hit_in = "True";
	  (_planes[pl]->GetHits()[iht])->set_name(candHit, hit_in);
	  
	  //pattern recognition chi2 
	  if ((_planes[pl]->GetHits()[iht] )->get_mu_prop() > 0.8 )
	    _recChi[0] = TMath::Max(Chi2[iht], _recChi[0]);
	  else
	    _recChi[1] = TMath::Min(Chi2[iht], _recChi[1]);
	  
	  
	  if ( nConsecHole > _recChi[2] )
	    _recChi[2] = nConsecHole;
	  
	  nConsecHole = 0;
	  
	} else nConsecHole++;
	
      } else {
	
	///assign as hadron
	    const dict::Key hadHit = "inhad";
	    const dict::Key had_in = "True";
	    (_planes[pl]->GetHits()[iht])->set_name(hadHit, had_in);
	    hads.push_back(_planes[pl]->GetHits()[iht]);
	    
      }
    }
    
  }//plane
  
  
  if ( nConsecHole > _max_consec_missed_planes ) {
    if ( _endProj ) check_forwards( seed, hits, muontraj );
    if ( muontraj.size() < _min_plane_prop*(double)(_nplanes-_badplanes) ){
      _failType = 6;
      return false;
    } else {
      sort_hits( hits, muontraj, hads );
      return true;
    }
  }
 
  if ( _endProj ) check_forwards( seed, hits, muontraj );
  
  _m.message("At the end of muon_extraction size in traj =",muontraj.size(), bhep::VERBOSE);
  return true;
}







//***********************************************************************
void event_classif::check_forwards(const State& seed, vector<cluster*>& hits,
				   Trajectory& muontraj){
  //***********************************************************************

  //Those not in 'blob' class get checked if there are any more hits
  //that can be added to the trajectory.
  

  _m.message("********Inside event_classif::check_forwards ", bhep::DETAILED);
 

  bool ok;
  int counter = 0;
  long ChiMin;
  double chi2[2];
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";
  


  //----------------------------------------------------
  for(int pl = _endLongPlane + 1; pl< _nplanes; pl++ ){
    
    if(_planes[pl]->GetNHits()!=2) break;
    
    
    for (int j=0;j<2;j++){
      ok = man().matching_svc().match_trajectory_measurement(muontraj, (*(_planes[pl]->GetHits()[j])), chi2[j]);
    }
    
    ChiMin = TMath::LocMin( 2, chi2);
      
    if ( ChiMin == 0 ){
     
      ok = man().fitting_svc().filter(*(_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -2]), seed, muontraj);
      
      if ( !ok ) (_planes[pl]->GetHits()[_planes[pl]->GetNHits()-2])->set_name(skipped,hit_in);
      else (_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -2])->set_name(candHit, hit_in);
      
      (_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -1])->set_name(skipped,hit_in);

    } else {
      ok = man().fitting_svc().filter(*(_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -1]), seed, muontraj);
      
      if ( !ok ) (_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -1])->set_name(skipped,hit_in);
      else (_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -1])->set_name(candHit, hit_in);

      (_planes[pl]->GetHits()[(_planes[pl]->GetNHits()) -2])->set_name(skipped,hit_in);
    }
  

    counter++;
  }
  
  
  _badplanes = counter;
  //-----------------------------------------------------  
}




//***********************************************************************
bool event_classif::invoke_cell_auto(vector<cluster*>& hits,
				     Trajectory& muontraj, vector<cluster*>& hads){
  //***********************************************************************
  //uses cellular automaton to try and retrieve more complicated events.
 
  _m.message("+++ Performing Cellular Automaton analysis +++",bhep::VERBOSE);
  
  bool ok;
  
  /// For radail search
  ///  if ( (_nplanes-_exclPlanes < _min_hits) && (_vRadCount.size() - _exclPlanes < _min_hits) ) return false;
  

  if (_nplanes-_exclPlanes < _min_hits)  return false;
  _m.message("_nplanes =",_nplanes,"   _exclPlanes=",_exclPlanes,"   _min_hits=",_min_hits, bhep::DETAILED);

  //TEST!!
  measurement_vector hit_meas;
  get_cluster_meas( hits, hit_meas );
  
  
  //vector of trajectories from CA
  std::vector<Trajectory*> trajs;
  ok = man().matching_svc().find_trajectories( hit_meas, trajs );
  _m.message(" in CA :: traj =",trajs.size(), bhep::DETAILED);



  ///if no track found
  if ( !ok || trajs.size() == 0) return false;
  
 
  if ( trajs.size() == 1 ) {
     if ((int)trajs[0]->size() >= _min_hits ){
 
      muontraj.reset();
      
      muontraj = *trajs[0];
      

      ///to assign candHit & hadrons from muontraj
      sort_hits( hits, muontraj, hads);
      
      ok = true;
    } else ok = false;
    
  } else {
    
    ok = sort_trajs( muontraj, trajs);
 
    if ( ok )
      sort_hits( hits, muontraj, hads);

  }
  
  
  ///delete the trajctories except the best one 
  if ( ok )
    delete_bad_trajs( muontraj, trajs );
  else  stc_tools::destroy( trajs );
  
  
  if (ok) _m.message("in CA:: traj info +++  , nmeas=",muontraj.size(), bhep::DETAILED);
  
  _m.message("+++ End of Cellular automaton +++",bhep::VERBOSE);
  
  return ok;
}





//***********************************************************************
void event_classif::get_cluster_meas(const vector<cluster*>& hits,
				     measurement_vector& meas)
  //***********************************************************************
{
  
  std::vector<cluster*>::const_iterator cIt;
  for (cIt = hits.begin()+_vertGuess;cIt != hits.end();cIt++)
    meas.push_back( (*cIt) );

}





//***********************************************************************
void event_classif::sort_hits(vector<cluster*>& hits,
			      Trajectory& muontraj, vector<cluster*>& hads){
  //***********************************************************************
  _m.message("*******event_classif::sort_hits*******", bhep::DETAILED);  

  bool inTraj;
  vector<Node*>::iterator inTrIt;
  inTrIt = muontraj.nodes().begin();

  const dict::Key candHit = "inMu";
  const dict::Key not_in = "False";
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";


  
  //Protection against double entries in hads in max_consec_planes fail case.
  if ( _intType != 5 ) hads.clear();

  
  //----------------------------------------
  for (int iht = 0; iht <(int)hits.size(); iht++){
    
    inTraj = false;
    if ( hits[iht]->names().has_key(candHit) )
      hits[iht]->set_name(candHit, not_in);
    
   
    ///check all the nodes of the muontraj
    for (int nd = 0; nd < (int)muontraj.nodes().size(); nd++){
    
      if ( hits[iht]->position() == muontraj.nodes()[nd]->measurement().position() ) {
	_m.message("I am here = ",hits[iht]->position()[2], bhep::DETAILED);   
	hits[iht]->set_name(candHit, hit_in);
	
	inTraj = true;
      }
     
    }
   
    ///if hits are not in Trajectory then assigned as hadron
    if ( !inTraj && !(hits[iht])->names().has_key(skipped) ){
      const dict::Key hadHit = "inhad";
      const dict::Key had_in = "True";
      hits[iht]->set_name(hadHit, had_in);
      hads.push_back(hits[iht] );
    }
    
  }
 

}




//***********************************************************************
void event_classif::delete_bad_trajs(const Trajectory& muontraj,
				     vector<Trajectory*>& trajs){
  //***********************************************************************

  vector<Trajectory*>::iterator dIt = trajs.begin();
  
  bool found = false;
  const dict::Key traj_index = "traj_no";
  
  while ( !found && dIt != trajs.end() ){
    if ( muontraj.quality( traj_index ) == (*dIt)->quality( traj_index ) ){
      found = true;
      
      trajs.erase( dIt );
    }
    dIt++;
  }
  stc_tools::destroy( trajs ); 
  

  /* /// At the end of the event destroy trajectories only valid for SINGLE TRACK FINDING CODE
  vector<cluster*> hits;
  while ( !found && dIt != trajs.end() ){
    if ( muontraj.quality( traj_index ) == (*dIt)->quality( traj_index ) ){
      found = true;

      trajs.erase( dIt );
    }
    dIt++;
  }
  if( _hitIt == hits.end()-1) stc_tools::destroy( trajs ); */ 
 

}




//***********************************************************************
bool event_classif::sort_trajs(Trajectory& muontraj, vector<Trajectory*>& trajs){
  //***********************************************************************
  //Reject possible trajectories based on number of hits, chi2, hits in common.
  
  bool ok;
  std::vector<Trajectory*> trajs2;
  
  ok = reject_small( trajs, trajs2);
  
  if ( ok ){
    if ( trajs2.size() == 1 ){
      
      muontraj.reset();
      
      muontraj = *trajs2[0];
      
      return ok;

    } else {

      std::vector<Trajectory*> trajs3;
      std::vector<double> Chis;
      
      ok = reject_high( trajs2, trajs3);
      
      //!!!!
      trajs2.clear();

      if ( ok ){
	if ( trajs3.size() == 1 ){
	  
	  muontraj.reset();
      
	  muontraj = *trajs3[0];
      
	  return ok;

	} else {
	  
	  ok = reject_final( trajs3, muontraj);
	  
	}

      }

    }

  }

  return ok;
}






//***********************************************************************
bool event_classif::reject_small(vector<Trajectory*>& trajs,
				 vector<Trajectory*>& trajs2){
  //***********************************************************************
  
  vector<Trajectory*>::iterator it1;

  for (it1 = trajs.begin();it1 != trajs.end();it1++){

    if ( (int)(*it1)->size() >= _min_hits )
      trajs2.push_back( (*it1) );

  }

  /* if ( trajs2.size() == 0 ) return false;

  return true;*/

  /// if trajs2 exist then only true
  if ( trajs2.size() == 0 ) return false;
  else  return true;
}





//***********************************************************************
bool event_classif::reject_high(vector<Trajectory*>& trajs,
				vector<Trajectory*>& trajs2){
  //***********************************************************************
  
  bool ok1, ok2 = false;
  double traj_no;
  const dict::Key traj_index = "traj_no";

  vector<Trajectory*>::iterator it1;
  State temp_seed;
  vector<cluster*> dummy_hits;
  double momErr;
  
  for (it1 = trajs.begin();it1 != trajs.end();it1++){
    
    Trajectory& temp_traj = *(*it1);
    traj_no = temp_traj.quality( traj_index );
    
    temp_traj.sort_nodes(RP::z, -1 );
    
    ok1 = get_patternRec_seed( temp_seed, temp_traj, dummy_hits);
    temp_traj.set_quality( traj_index, traj_no );
    
    if ( ok1 && temp_traj.quality() < _chi2_max ){

      (*it1)->set_quality( temp_traj.quality() );
      
      const dict::Key relErr = "relErr";
      momErr = (double)(temp_seed.matrix()[5][5] / temp_seed.vector()[5]);
      (*it1)->set_quality( relErr, momErr);

      trajs2.push_back( (*it1) );

      ok2 = true;
    }
    
  }

  sort( trajs2.begin(), trajs2.end(), chiSorter() );

  return ok2;
}






//***********************************************************************
bool event_classif::reject_final(vector<Trajectory*>& trajs,
				 Trajectory& muontraj){
  //***********************************************************************
  //Accept lowest local chi and then any others with few coincidences
  //with this trajectory then choose the longest/lowest error track.
  
  std::vector<Trajectory*> temp_trajs;
  temp_trajs.push_back( trajs[0] );

  vector<Trajectory*>::iterator it1, it2;

  double coincidence = 0;

  for (it1 = trajs.begin()+1;it1 != trajs.end();it1++){
    
    vector<Node*>& node2 = (*it1)->nodes();
    it2 = temp_trajs.begin();

    while ( it2 != temp_trajs.end() && coincidence < _max_coincedence ){

      vector<Node*>& node1 = (*it2)->nodes();

      coincidence = compare_nodes( node1, node2);

      it2++;
    }
    
    if ( coincidence < _max_coincedence )
      temp_trajs.push_back( (*it1) );

  }
  trajs.clear();
  
  if ( temp_trajs.size() == 1 ) {
    
    muontraj.reset();

    muontraj = *temp_trajs[0];

  } else {

    /// sort trajectory by hits instead of length
     sort( temp_trajs.begin(), temp_trajs.end(), sortTrajByHits() );
     muontraj = *temp_trajs[0];
    // select_trajectory( temp_trajs, muontraj);

  }
  temp_trajs.clear();

  return true;
}







//***********************************************************************
double event_classif::compare_nodes(const vector<Node*>& n1,
				    const vector<Node*>& n2){
  //***********************************************************************
  
  std::vector<Node*>::const_iterator nIt1, nIt2;
  double counter = 0;

  for (nIt1 = n1.begin();nIt1 != n1.end();nIt1++){
    for (nIt2 = n2.begin();nIt2 != n2.end();nIt2++){

      if ( (*nIt1)->measurement().position() == (*nIt2)->measurement().position() )
	counter++;

    }
  }

  return counter / (double)n2.size();
}

//***********************************************************************
void event_classif::select_trajectory(vector<Trajectory*>& trajs,
				      Trajectory& muontraj){
  //***********************************************************************
  
  int length[(const int)trajs.size()];
  double Err[(const int)trajs.size()];
  
  muontraj.reset();
  
  const dict::Key relErr = "relErr";
  
  for (int it1 = 0;it1 <(int)trajs.size();it1++){

    length[it1] = (int)trajs[it1]->length();
    Err[it1] = trajs[it1]->quality( relErr );
    
  }

  int longest = (int)TMath::LocMax( (int)trajs.size(), length);
  int safest = (int)TMath::LocMin( (int)trajs.size(), Err);
  
  if ( longest == safest )
    muontraj = *trajs[longest];
  else {

    if ( length[longest] - length[safest] > 5)
      muontraj = *trajs[longest];
    else muontraj = *trajs[safest];

  }

}

//***********************************************************************
void event_classif::set_branches(){
  //***********************************************************************
  
  _likeTree->Branch("truInt",&_truInt,"Int/I");
  _likeTree->Branch("nhits", &_nhit, "nhits/I");
  _likeTree->Branch("VisEng", &_visEng, "visEng/D");
  _likeTree->Branch("nplanes", &_nplanes, "nplanes/I");
  _likeTree->Branch("freePlanes", &_freeplanes, "freeplanes/I");
  _likeTree->Branch("occupancy", &_occ, "occ[nplanes]/I");
  _likeTree->Branch("planeEnergy", &_plEng, "plEng[nplanes]/D");
  _likeTree->Branch("nhitTraj", &_trajhit, "nhitTraj/I");
  _likeTree->Branch("TrajPur", &_trajpur, "trajpur/D");
  _likeTree->Branch("EngTraj", &_trajEng, "engTraj/D");
  _likeTree->Branch("EngTrajPlane", &_trajEngPlan, "engTrajPlane[nplanes]/D");
  _likeTree->Branch("HitsInTrajClust",&_trclusthits,"nhittrajclust[nplanes]/I");

}

//***********************************************************************
void event_classif::output_liklihood_info(const vector<cluster*>& hits){
  //***********************************************************************
  _m.message("I am in output_liklihood_info",bhep::DETAILED);
  
  bool multFound = false;
  
  _nhit = 0;
  _freeplanes = 0;
  _visEng = 0;
  _trajhit = 0;
  _trajpur = 0;
  _trajEng = 0;

  for (int ipl = 0;ipl<_nplanes;ipl++){
    _trajEngPlan[ipl] = 0;
    _trclusthits[ipl] = 0; 
  }
  
  
  
  for (int pl = _nplanes -1; pl>=0; pl--){
    
    _nhit += _planes[pl]->GetNHits();
    _visEng += _planes[pl]->GetEng();
    
    _occ[pl] = _planes[pl]->GetNHits();
    _plEng[pl] = _planes[pl]->GetEng();
    
    if ( !multFound && _planes[pl]->GetNHits() == 1 )
      _freeplanes++;
    else multFound = true;

    _m.message("in outlke:: _nhit=",_nhit," _visEng=",_visEng,"  _freeplanes=",_freeplanes,"   _occ=",_occ[pl],"  _plEng=",_plEng[pl], bhep::DETAILED);   
  }
  //---------------------------------------------------------------
  
}


//***********************************************************************
void event_classif::traj_like(const vector<cluster*>& hits){
  //***********************************************************************
  
  const dict::Key Edep = "E_dep";
  const dict::Key candHit = "inMu";
  vector<cluster*>::const_iterator hitIt3;
  //vector<Node*>::const_iterator trIt1 = muontraj.nodes().begin();
  vector<double>::iterator planIt;

  int counter = _nplanes - 1;
  
  for (hitIt3 = hits.begin();hitIt3 != hits.end();hitIt3++){

    ///if it is a candidate muon
    if ( (*hitIt3)->names().has_key( candHit ) ){
      if ( (*hitIt3)->name( candHit ).compare("True") == 0 ){
	_trajhit++;
	//_trajEng += bhep::double_from_string( (*hitIt3)->name( Edep ) ) * GeV;
	_trajEng += (*hitIt3)->get_eng() * MeV;
	_trajEngPlan[counter] = (*hitIt3)->get_eng() * MeV;
	// if ( (*hitIt3)->name("MotherParticle").compare("mu+") == 0
	// 	     || (*hitIt3)->name("MotherParticle").compare("mu-") == 0 )

	///if the muon candidate is true muon
	if ( (*hitIt3)->get_mu_prop() > 0.8 )//still to be decided.
	  _trajpur++;
	_trclusthits[counter] = (*hitIt3)->get_nVox();
	_m.message("_trajhit=",_trajhit,"  _trajEng=",_trajEng,"  counter=",counter," _trajEngPlan=",_trajEngPlan[counter],"    _trclusthits= ",_trclusthits[counter], bhep::DETAILED);
	counter--;
      }
    }
  }
  _trajpur /= _trajhit;

 
 
}

//***********************************************************************
void event_classif::out_like(){
  //***********************************************************************
 
  /// int fillcatch = _likeTree->Fill();
  _likeTree->Fill();
  
}




//***********************************************************************
//Temp function for likelihoods with tru interaction in tree.
void event_classif::set_int_type(const string name){
  //***********************************************************************

  if ( name=="CCQE" )
    _truInt = 1;
  else if ( name=="NCQE" )
    _truInt = 2;
  else if ( name=="CCDIS" )
    _truInt = 3;
  else if ( name=="NCDIS" )
    _truInt = 4;
  else if ( name=="1piRes" )
    _truInt = 5;
  else if ( name=="miscRes" )
    _truInt = 6;
  else if ( name=="eEl" )
    _truInt = 7;
  else if ( name=="muINVe" )
    _truInt = 7;
  else if ( name=="unknown" )
    _truInt = -1;
  else
    _truInt = 8;

}


//***********************************************************************
double event_classif::correctEdep(double edep, double X, double Y, double Z)
  //***********************************************************************
{

  double sum1 = 0;
  
  double slope = _octGeom==1 ? (_detY - _detX*tan(atan(1)/2.))/
    (_detY*tan(atan(1)/2.) - _detX) : -1.0;
  if(_octGeom==0 && Z > (_vdetZ - _detZ)/2.){ // Assume a cylindrical geometry
    double xedge = _detX * sqrt(1 - 4.*pow(Y/_detY, 2))/2.;
    double yedge = _detY * sqrt(1 - 4.*pow(X/_detX, 2))/2.;
    sum1 =  exp( -(xedge - fabs(X))/_WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
  }
  else if(_octGeom==0 && Z < (_vdetZ - _detZ)/2.){ // Assume a cylindrical geometry
    double xedge = _detX/2.;
    double yedge = _detY/2.;
    sum1 =  exp( -(xedge - fabs(X))/_WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
  }
  else if(_octGeom==1) {
    //need to take into account drift distance to closest and furthest edge.
    if((fabs(X) < _detY*tan(atan(1)/2.)/2 &&
	fabs(Y) < _detX*tan(atan(1)/2.)/2 ) ){
      sum1 =  exp( -(_detX/2 - fabs(X))/_WLSAtten );
      sum1 += exp( -(_detX/2 + fabs(X))/_WLSAtten );
      sum1 += exp( -(_detY/2 - fabs(Y))/_WLSAtten );
      sum1 += exp( -(_detY/2 + fabs(Y))/_WLSAtten );
    }
    else if(fabs(X) > _detY*tan(atan(1)/2.)/2 &&
	    fabs(Y) < _detX*tan(atan(1)/2.)/2  ){
      double xedge = _detX/2 + (fabs(Y) - 
				_detX/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge  - fabs(X))/_WLSAtten );
      sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
      sum1 += exp( -(_detY/2 - fabs(Y))/_WLSAtten );
      sum1 += exp( -(_detY/2 + fabs(Y))/_WLSAtten );
    }
    else if(fabs(X) < _detY*tan(atan(1)/2.)/2 &&
	    fabs(Y) > _detX*tan(atan(1)/2.)/2  ){
      double yedge = _detY/2 + (fabs(X) - 
				_detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(_detX/2 - fabs(X))/_WLSAtten );
      sum1 += exp( -(_detX/2 + fabs(X))/_WLSAtten );
      sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
      sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
    }
    else if(fabs(X) > _detY*tan(atan(1)/2.)/2 &&
	    fabs(Y) > _detX*tan(atan(1)/2.)/2  ){
      double xedge = _detX/2 + (fabs(Y) - 
				_detX/2.*tan(atan(1)/2.))*slope;
      double yedge = _detY/2 + (fabs(X) - 
				_detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge - fabs(X))/_WLSAtten );
      sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
      sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
      sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
    }
  }
  else if(_octGeom==2) {
    //need to take into account drift distance to closest and furthest edge.
    if((fabs(X) < _detY*tan(atan(1)/2.)/2. + _detX/4.&&
	fabs(Y) < _detY*tan(atan(1)/2.)/2. ) ){
      sum1 =  exp( -(_detX/2 - fabs(X))/_WLSAtten );
      sum1 += exp( -(_detX/2 + fabs(X))/_WLSAtten );
      sum1 += exp( -(_detY/2 - fabs(Y))/_WLSAtten );
      sum1 += exp( -(_detY/2 + fabs(Y))/_WLSAtten );
    }
    else if(fabs(X) > _detY*tan(atan(1)/2.)/2. + _detX/4. &&
	    fabs(Y) < _detY*tan(atan(1)/2.)/2.  ){
      double xedge = _detX/2 + (fabs(Y) - 
				_detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge  - fabs(X))/_WLSAtten );
      sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
      sum1 += exp( -(_detY/2 - fabs(Y))/_WLSAtten );
      sum1 += exp( -(_detY/2 + fabs(Y))/_WLSAtten );
    }
    else if(fabs(X) < _detY*tan(atan(1)/2.)/2 + _detX/4. &&
	    fabs(Y) > _detY*tan(atan(1)/2.)/2  ){
      double yedge = _detX/2 + (fabs(X) - 
				_detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(_detX/2 - fabs(X))/_WLSAtten );
      sum1 += exp( -(_detX/2 + fabs(X))/_WLSAtten );
      sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
      sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
    }
    else if(fabs(X) > _detY*tan(atan(1)/2.)/2 + _detX/4. &&
	    fabs(Y) > _detY*tan(atan(1)/2.)/2  ){
      double xedge = _detX/2 + (fabs(Y) - 
				_detX/2.*tan(atan(1)/2.))*slope;
      double yedge = _detX/2 + (fabs(X) - 
				_detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge - fabs(X))/_WLSAtten );
      sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
      sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
      sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
    }
  }
  else if(_octGeom == -2 || _octGeom == 3){ // Assume a square cross-section
    double xedge = _detY/2.;
    double yedge = _detY/2.;
    sum1 =  exp( -(xedge - fabs(X))/_WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
  }
  double corrEng = 4*edep/sum1;
  return corrEng;
}

/***************************************************************************************/
double event_classif::MomentumFromDeflection(const Trajectory& traj, int firsthit){
  /***************************************************************************************/

  // Calculate the momentum from the deflection angle between hits.
  // The momentum is proportional to the magnetic field * change in angle.

  double fac = 1;//Should be prop to q*b/m measure emperically.

  EVector Z = EVector(3,0); Z[2] = 1;

  double sumRelP =0;

  //double num_nodes = traj.size() -3 - firsthit;

  double num_nodes = traj.size() -1;

  cout<<num_nodes<<endl;

  // for (int ipoint=firsthit;ipoint <=num_nodes;ipoint++)
  //for (int ipoint=num_nodes;ipoint >= 2;ipoint--)
  // for (int ipoint=num_nodes-1;ipoint >= num_nodes-4;ipoint--)
  //{
      //EVector pos0 = traj.node(ipoint).measurement().position();
      //EVector pos1 = traj.node(ipoint + 1).measurement().position();
      //EVector pos2 = traj.node(ipoint + 2).measurement().position();

  int ipoint = num_nodes -1;

      EVector pos0 = traj.node(ipoint).measurement().position();
      EVector pos1 = traj.node(ipoint - 1).measurement().position();
      EVector pos2 = traj.node(ipoint - 2).measurement().position();

      EVector B0 = _geom.getBField(pos0);
      EVector B1 = _geom.getBField(pos1);

      cout<<pos0[2]<<" "
	  <<pos1[2]<<" "
	  <<pos2[2]<<" "
	  <<B0[0]<<" "
	  <<B1[0]<<endl;
      /*
      if(dot(B0,B1) < 0)
	{
	  //sumRelP = 0;
	  //num_nodes = ipoint -firsthit;
	  num_nodes = num_nodes -ipoint;
	  break;
	}
      */
      // Transverse magnetic field component
      EVector bT = crossprod(Z, B0);
      double bTcomp = sqrt(dot(bT,bT));


      EVector dPos0 = pos1-pos0;
      dPos0 /= sqrt(dot(dPos0,dPos0));
      EVector dPos1 = pos2-pos1;
      dPos1 /= sqrt(dot(dPos1,dPos1));

      double deflection_angle =sqrt(1-dot(dPos0,dPos1));

      cout<<"ang: "<< acos(dot(dPos1,Z)/sqrt(dot(dPos1,dPos1)))<<" "
	  <<acos(dot(dPos0,Z)/sqrt(dot(dPos0,dPos0)))<<endl;
      
      cout<<"diff ang: "<< acos(dot(dPos1,Z)/sqrt(dot(dPos1,dPos1))) -
	acos(dot(dPos0,Z)/sqrt(dot(dPos0,dPos0)))<<endl;

      cout<<"angle: "<<deflection_angle<<endl;

      sumRelP += fac * deflection_angle;// * bTcomp;
      //  }

  //double p = sumRelP/(num_nodes);
  double p = sumRelP/4;

  cout<<"momentum guess from deflection: "<<p<<endl;

  return p;
}


/***************************************************************************************/
double event_classif::RangeMomentum(double length, double nodeZ){
  /***************************************************************************************/
  std::map<dict::Key,vector<double> > moduleDataMap = _geom.getModuleDataMap();
  double p = 0;

  for (std::map<dict::Key,vector<double> >::iterator it=moduleDataMap.begin();
       it!=moduleDataMap.end(); ++it)
    {
      double module_pos = it->second[0];
      double module_half_size = it->second[1];
      double wFe = it->second[2];

      //std::cout<<"Fitter "<<module_pos<<" "<<module_half_size<<" "
      //     <<wFe<<" "<<nodeZ<<" "<<length<<std::endl;

      std::cout<<"Fitter "<<nodeZ<<" "<<length<<" "<<wFe<<" "<<module_pos<<std::endl;

      //Sanity checking not outside range forward 
      if(!((module_pos-module_half_size)<(nodeZ+length)))
	{
	  continue;
	  std::cout<<"Fitter not in range"<<std::endl;
	  
	}

      // Through the whole module
      else if((nodeZ + length) > (module_pos + module_half_size))
	{
	   p += 2*module_half_size * wFe;
	   std::cout<<"p+="<<2*module_half_size * wFe<<std::endl;
	   std::cout<<"Fitter through"<<std::endl;
	}
      // Stop in the module
      else if((nodeZ + length) < (module_pos + module_half_size))
	{
	  p+=((length + nodeZ) - (module_pos - module_half_size))*wFe;
	  //p+=(module_half_size + module_pos -(length+nodeZ))*wFe;
	  std::cout<<"p+="<<((length + nodeZ) - (module_pos - module_half_size))*wFe<<std::endl;
	  std::cout<<"Fitter stop"<<std::endl;
	}	 
    }

  //std::cout<<"In Map p "<<p<<std::endl;  


  // Computing the momentum of a track based on the range
  // Parameters calculated from "Atomic Data and Nuclear Data Tables 78, 183-356 (2001)"
  //double p = (_FeWeight*(0.011844*GeV * pow((length/cm),1.03235))
  //      + (1- _FeWeight)*(0.0023705067*GeV* pow(length/cm,1.00711577)));
  return p;
  // Should be compared to G4 simulation directly.
}




//***********************************************************************
void event_classif::fill_traj_info( Trajectory& muontraj) {
  //***********************************************************************
  _m.message("++++event_classif::fill_traj_info ", bhep::VERBOSE);
  
  muontraj.set_quality("intType",_intType);
  muontraj.set_quality("failType",_failType);
  muontraj.set_quality("nplanes",_nplanes);
  muontraj.set_quality("freeplanes",_freeplanes);
  muontraj.set_quality("badplanes",_badplanes);
  muontraj.set_quality("longestSingle",_longestSingle);
  muontraj.set_quality("radialLongest",_radialLongest);
  muontraj.set_quality("xtent",_Xtent);
  muontraj.set_quality("vertZ",_vertGuessZ);
  muontraj.set_quality("lowPt",_lowPt);
  
  /// for CA tracks it is not sure the free section belongs to the same traj or not
  if((int)muontraj.quality("intType")!=5) muontraj.set_quality("lastIso", _lastIso);
  else muontraj.set_quality("lastIso", 0);
  
  _m.message(" event_clss failType =",_failType, bhep::VERBOSE);
}

 
//***********************************************************************
void event_classif::reset_traj_info() {
  //***********************************************************************
  _m.message("++++event_classif::reset_traj_info ", bhep::VERBOSE);
  //reset the containers
  _nplanes = 0;
  _meanOcc = 0;
  _badplanes = 0;
  _longestSingle = 0;
  _endProj = false;
  _vRadCount.clear();
  _radialLongest = 0;
  _freeplanes = 0;
  _endLongSing = 0;
  _endLongPlane = 0;
  _lastIso = 0;
  _intType = 0;
  _Xtent = 0;
  _vertGuess = 0;
  _vertGuessZ = 0;
  _failType = 0;

  _nRadPlanes=0;
  _planeEnd.reset();

  _planes.clear();
  _radPlanes.clear();

  stc_tools::destroy(_planes);
  stc_tools::destroy(_radPlanes);
  

}
/// Radial search not yet incorporatedreset

//***********************************************************************
bool event_classif::get_radial_occupancy(vector<cluster*>& hits)
 ///***********************************************************************
{
   _m.message("++++ Calculating Radial occupancies ++++",bhep::VERBOSE);
 
  
  bool ok = true;
  
 
  int single_count = 0, planeIndex=0;
  double  testX, curX, testY, curY, testZ, diffR, testR; /// EngPlane = 0, 
  size_t hits_used = 0, imeas = 0;
  _radialLongest = 0;
  _radialFree = 0;
  
  vector<cluster*> hits_temp = hits;
  size_t nHits = hits_temp.size();
  
  
  ///sort hits by R independent of Z  
  // sort(hits_temp.begin(), hits_temp.end(), sortHitsByR());
  // cout<<"  radial_hits end Z="<<hits_temp.back()->position()[2]<<endl;  


  ///start count radial occupancy
  do {
       
    testX = hits_temp[imeas]->position()[0];
    testY = hits_temp[imeas]->position()[1];
    testZ = hits_temp[imeas]->position()[2];
    testR = sqrt(testX*testX + testY*testY);
    hits_used++;

    ///create plane info
    plane_info* radPlane = new plane_info(planeIndex, testZ, testR, _infoStore);
    radPlane->AddHit(hits[imeas]);
    
    ///calculate the z position which is the current z for hits 1 -> total no of hits in the cluster 
    for (size_t i = hits_used;i <nHits;i++) {
      curX = hits_temp[i]->position()[0];
      curY = hits_temp[i]->position()[1];
      diffR = sqrt(pow(testX - curX, 2) + pow(testY - curY, 2));
      ///tolerance 5mm in R
      if (diffR <= sqrt(2.)*_voxEdge) {
	_m.message("testZ = ",testZ,"   diffR = ",diffR,bhep::VERBOSE);
	testX = hits_temp[i]->position()[0];
	testY = hits_temp[i]->position()[1];
	testR = sqrt(testX*testX + testY*testY);
	//add the hit
	radPlane->AddHit(hits_temp[i]);

	hits_used++;
      } else break;
	
    }
      
    //_vRadCount.push_back(count);
    _radPlanes.push_back(radPlane);
      
 
    ///if single hit plane then increase the single_count, otherwise assign it 0
    if ( radPlane->GetNHits() == 1 )
      single_count++;
    else {
      if ( single_count >= _radialLongest ){
	_radialLongest = single_count;
	_endLongSing = imeas -1;
	_endLongPlane = _radPlanes.size() -2;
	
      }
      single_count = 0;
    }
    if ( imeas == nHits -1 && single_count != 0 && radPlane->GetNHits() == 1 ){
      if ( single_count >= _radialLongest ){
	_radialLongest = single_count;
	_endLongSing = imeas;
	_endLongPlane = _radPlanes.size() -1;

	
      }
    }
    
    //
    planeIndex++;
    ///increase the measurement count and also mean0cc
    imeas += radPlane->GetNHits();
   

    ///occupancy
    //totalOcc += radPlane->GetNHits();  
    
  } while (hits_used != nHits);

  hits_temp.clear();
  
  ///total no of planes
  _nRadPlanes = (int)_radPlanes.size();
  
  
  /// Mean occupancy
  if ( _nRadPlanes == 0 ) return false;
  // _meanOcc = totalOcc/ (double)_nRadPlanes;
  
  return ok;
}

//***********************************************************************
bool event_classif::LowMomentumExtraction(vector<cluster*>& hits,
				    Trajectory& muontraj, vector<cluster*>& hads) {
  //***********************************************************************

  _m.message("++++ event_classif::low_momentum_extraction +++++++++++++", bhep::VERBOSE);
 
  //Calculate the momentum and the charge for short tracks.

  bool ok =false;
  int nMeas = muontraj.size();
  double firstNodeZ = muontraj.nodes()[nMeas-1]->measurement().position()[2];
  double pathlength=muontraj.nodes()[0]->measurement().position()[2] - firstNodeZ;

  double p = RangeMomentum(pathlength,firstNodeZ);
  
  double q = CalculateCharge(muontraj);

  // Give the momentum a charge.
  p *= q;

  cout<<"For low momentum: "<<p<<" "<<q<<endl;

  for(unsigned int cnt = 0; cnt < muontraj.nodes().size(); cnt++)
    {
      muontraj.nodes()[cnt]->set_status(RP::fitted);
      
    }
  
  // How give this to the track without using the fitter? Can we give it to the fitter given so few hits?
  // How pass through fitter?

  //State& seedState = muontraj.nodes()[0]->state();
  State seedState;
  //EVector v = seedState.vector();
  EVector V(6,0);
  EMatrix M(6,6,0);

  V[5] = 1/p;
  //EMatrix C0 = seedState.matrix();
  seedState.set_name(RP::particle_helix);
  seedState.set_name(RP::representation,RP::slopes_curv_z);
  seedState.set_hv(RP::sense,HyperVector(V,M,RP::x));
  seedState.set_hv(HyperVector(V,M,RP::slopes_curv_z));
  
  //State currstate = traj.state(traj.first_fitted_node());
  //muontraj.nodes()[0]->set_state(seedState);
  muontraj.nodes()[muontraj.first_fitted_node()]->set_state(seedState);
  muontraj.set_quality("fitted",true);
  muontraj.set_quality("initialqP",p);
  muontraj.set_quality("fitcheck", 4);
  //muontraj.set_quality("lowPt",1);

  _lowPt = 1;

  ok = true;

  /*
  State patternSeed;
  
  ok = get_patternRec_seed( patternSeed, muontraj, hits);
  if(!ok)  _m.message(" get_patternRec_seed not ok",bhep::DETAILED); 
  
  if ( ok )
    ok = perform_muon_extraction( patternSeed, hits, muontraj, hads);
  if(!ok) _m.message("perform_muon_extraction not ok",bhep::DETAILED);

  ///assign the seed state
  if ( ok )
    _seedState = patternSeed;
  
  if(ok) _m.message(" event_class: traj nmeas=",muontraj.size()," intType =",_intType,"  && PR seed is=",_seedState,bhep::DETAILED);

  /// to get all the PR seed for reseed_traj inside fitter
  
  if(ok) _vPR_seed.push_back(_seedState);
  
  */
  return ok;
}

//***********************************************************************
double event_classif::CalculateCharge(Trajectory& track) {
  //***********************************************************************

  int fitcatcher;
  int nMeas = track.size();

  if (nMeas > 4) nMeas = 4;

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];/// r[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;
  double pathlength=0;/// tminR=9999.9;

  double firstNodeZ = track.nodes()[nMeas-1]->measurement().position()[2];
  
  pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
  pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
  z[nMeas-1] = track.nodes()[nMeas-1]->measurement().position()[2];
  pos[2] = 0.0;

  for (int iMeas = nMeas-2;iMeas >= 0;iMeas--){
    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];
    // get the b-field from the previous step
    EVector B = _geom.getBField(pos);
    u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
    pos[0] = x[iMeas];  pos[1] = y[iMeas];   pos[2] = z[iMeas];
    EVector dR = EVector(3,0);
    dR[0] = x[iMeas+1] - x[iMeas];
    dR[1] = y[iMeas+1] - y[iMeas];
    dR[2] = z[iMeas+1] - z[iMeas];
    double dq = 0.;
    //pathlength += dR.norm();
    if(iMeas < nMeas-3){
      EVector dR1 = EVector(3,0);
      dR1[0] = x[iMeas+2] - x[iMeas+1];
      dR1[1] = y[iMeas+2] - y[iMeas+1];
      dR1[2] = z[iMeas+2] - z[iMeas+1];
      EVector ddR = dR1 + dR;
      dq = dot(ddR, crossprod(Z,B))/ (crossprod(Z,B).norm());
    }
    if(pdR != 0)
      if(// minR < R && dR/fabs(dR) == -pdR/fabs(pdR) && 
	 (x[iMeas]/fabs(x[iMeas]) != x[iMeas+1]/fabs(x[iMeas+1]) ||
	  y[iMeas]/fabs(y[iMeas]) != y[iMeas+1]/fabs(y[iMeas+1]))){ 
	// Sign change and minimum momentum
	minR = pos.norm();
	minindex = iMeas;
	sumdq = 0.0;
      }
    pdR = dR.norm();
    sumdq += dq;
  }

  TGraph *gr = new TGraph((const int)minindex, z, u);

  TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun->SetParameters(0.,0.001,0.001);

  fitcatcher = gr->Fit("parfit2", "QN");

  double qtilde = -0.3*pow(1+fun->GetParameter(1),3./2.)/
  (2*fun->GetParameter(2));
  int meansign = (int)(qtilde/fabs(qtilde));
  //delete fun;
  //return qtilde;
  //int meansign = sumdq/fabs(sumdq);
  
  return meansign;
}

/*
bool event_classif::get_patternRec_seed_low_track(State& seed, Trajectory& muontraj,
					vector<cluster*>& hits) {
 

  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);
  
  ///x,y,z values of the 1st node
  V[0] = muontraj.nodes()[0]->measurement().vector()[0];
  V[1] = muontraj.nodes()[0]->measurement().vector()[1];
  V[2] = muontraj.nodes()[0]->measurement().position()[2];

  //direction
  double dqtot = fit_parabola( V, muontraj);

  //Pos. gradient in X is ~negative curvature.
  // Assume that the R-Z plane is the bending plane for the toroidal field
  double VRad =  dqtot/fabs(dqtot);

  cout<<"In get_patternRec_seed: "<<V[5]<<" "<<1/V[5]<<endl;
  
  //Errors
  M[0][0] = M[1][1] = 15.*cm*cm;
  M[2][2] = EGeo::zero_cov()/2;
  M[3][3] = M[4][4] = 1.5;
  M[5][5] = pow(V[5],2)*4;

  //Sense
  V2[0] = 1;
  
  //Seedstate fit properties
  seed.set_name(RP::particle_helix);
  seed.set_name(RP::representation,RP::slopes_curv_z); 
  seed.set_hv(RP::sense,HyperVector(V2,M2,RP::x));
  seed.set_hv(HyperVector(V,M,RP::slopes_curv_z));
  
  // std::cout<<seed.hv().representation()<<std::endl;

  // std::cout<<"Is good here"<<std::endl;
  //  man().model_svc().model(RP::particle_helix).representation(RP::slopes_curv_z)
  //     .convert(seed,RP::default_rep);
  
  bool ok = perform_kalman_fit( seed, muontraj);
  
  if ( !ok )
    _failType = 5;
  
  return ok;
}
*/

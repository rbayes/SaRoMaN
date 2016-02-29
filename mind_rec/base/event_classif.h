/* -*- mode: c++ -*- */
#ifndef _event_classif___
#define _event_classif___

#include <mind/MINDsetup.h>
#include <mind/SetupSk.h>
#include <mind/MINDfitman.h>
#include <mind/cluster.h>

#include <mind/super_fit.h>

#include <recpack/RecPackManager.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/ParticleState.h>
#include <recpack/CellularAutomaton.h>

#include <bhep/gstore.h>

#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>


///
#include <mind/plane_info.h>

class event_classif : public super_fit{
  
public:
  
  event_classif();
  
  virtual ~event_classif();
  
  //-------------- main functions --------------//
  //void initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, Setup& det, double wFe);
  void Initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, double wFe, MINDsetup* geom);
  bool Execute(const vector<cluster*>& hits,
	       vector<Trajectory*>& vtrajs, vector<cluster*>& hads);
  void Finalize();
  //-------------------------------------------//
  
  //Grabbers for monitoring etc.
  int get_int_type(){ return _intType; }
  int get_vertex(){ return _vertGuess; }
  int get_fail_type(){ return _failType; }
  std::vector<EVector>& get_PatRec_Chis(){ return _recChi; }
  State& get_patRec_seed(){ return _seedState; }
  std::vector<State>& get_patRec_seed_vector(){ return _vPR_seed; }////
  int get_last_iso(){ return _lastIso; }
  double get_vis_eng(){ return _visEng; }
  int get_planes(){ return _nplanes; }
  int get_free_planes(){ return _freeplanes; }
  double get_xtent(){ return _Xtent; }
  

 
 
//Getters.    
  // Trajectory& get_best_traj(){ return best_traj; }
  
  bool get_plane_occupancy(vector<cluster*>& hits);
  bool get_radial_occupancy(vector<cluster*>& hits);///
  void assess_event(vector<cluster*>& hits);
 
  
  //Tempory for likelihoods.
  void set_int_type(const string name);
  //R
  double correctEdep(double edep, double X, double Y, double Z);

  //double RangeMomentum(double length,double nodeZ);
  //double MomentumFromCurvature(const Trajectory& traj, int firsthit = 0);
  
protected:
  
  void readParam();
  //void set_extract_properties(Setup& det);
  void reset();
  
  //Functions to be performed on CC mu candidates.  
  bool chargeCurrent_analysis(vector<cluster*>& hits,
			      Trajectory& muontraj, vector<cluster*>& hads);
  void fill_traj_info(Trajectory& muontraj);
  void reset_traj_info();
  //int exclude_backwards_particle();
  bool muon_extraction(vector<cluster*>& hits,
		       Trajectory& muontraj, vector<cluster*>& hads);
  
  // Handle low momentum tracks.
  bool LowMomentumExtraction (vector<cluster*>& hits,
		       Trajectory& muontraj, vector<cluster*>& hads);

  // Calculate the charge of a track using scattering angles.
  //double CalculateCharge(Trajectory& track);
  //double CalculateCharge2(Trajectory& track);


  bool get_patternRec_seed(State& seed, Trajectory& muontraj, vector<cluster*>& hits);
  double fit_parabola(EVector& vec, Trajectory& track);
  void set_de_dx(double mom);
  bool perform_kalman_fit(State& seed, Trajectory& track);
  bool perform_muon_extraction(const State& seed, vector<cluster*>& hits,
			       Trajectory& muontraj, vector<cluster*>& hads);
  void check_forwards(const State& seed, vector<cluster*>& hits,
		      Trajectory& muontraj);

  // void use_mini_cellAuto(const int occ, Trajectory& muontraj);
  //specific functions using cellular automaton.
  bool invoke_cell_auto(vector<cluster*>& hits,
			Trajectory& muontraj, vector<cluster*>& hads);
  void sort_hits(vector<cluster*>& hits, Trajectory& muontraj, vector<cluster*>& hads);
  void get_cluster_meas(const vector<cluster*>& hits, measurement_vector& meas);
  void delete_bad_trajs(const Trajectory& muontraj, vector<Trajectory*>& trajs);
  bool sort_trajs(Trajectory& muontraj, vector<Trajectory*>& trajs);
  bool reject_small(vector<Trajectory*>& trajs, vector<Trajectory*>& trajs2);
  bool reject_high(vector<Trajectory*>& trajs, vector<Trajectory*>& trajs2);
  bool reject_final(vector<Trajectory*>& trajs, Trajectory& muontraj);
  double compare_nodes(const vector<Node*>& n1, const vector<Node*>& n2);
  void select_trajectory(vector<Trajectory*>& trajs, Trajectory& muontraj);
  //


  RecPackManager& man(){
    return MINDfitman::instance().manager();}
  
  bhep::gstore _infoStore;
  
  bhep::messenger _m;
  
  MINDsetup _geom;

  int _lowPt;

  //Members to store plane occupancy and mean energy.
  int _nplanes;
  int _maxNtraj;
  double _meanOcc;
  vector<int> _hitsPerPlane;
  vector<double> _energyPerPlane;
  // vector<double> _planeZ;
  double _tolerance; //required 'closeness' to be considered in plane.
  double _voxEdge; //Voxel edge size to check for badly reconstructed points at end.

  std::vector<plane_info*> _planes;///
  std::vector<plane_info*> _radPlanes;///
  
  //integer for type candidate (NC etc.)
  int _intType;
  int _lastIso;
  plane_info _planeEnd;
 
  
  int _vertGuess;
  double _vertGuessZ;
  int _exclPlanes;
  int _badplanes;
  int _longestSingle;//length (in hits) of longest 'free' section.
  int _endLongSing;//End point of the above (position in hit vector).
  int _endLongPlane;//Plane position of above;
  double _maxBlobSkip;//max proportion of planes for basic skip.
  double _minBlobOcc;
  bool _endProj;//bit to tell if a forwards projection is needed
  
  //Monitoring variables.
  int _failType;
  std::vector<EVector> _recChi;
  State _seedState;

  double _FeWeight;
  double _detLength
;
  int _max_consec_missed_planes;
  int _min_seed_hits;
  int _min_check;
  int _min_hits;
  double _min_plane_prop;
  double _chi2_max;
  double _max_coincedence;
  double _pieceLength;
  
  bool _isTASD;
  //
  int _max_hits;
  int _max_nmult;

  //Output Likilihood info?
  bool _outLike;
  TFile *_outFileEv;
  TTree *_likeTree;

  /// for radial search
  vector<int> _vRadCount ;
  int _radialLongest;///
  int _radialFree;///
  int _nRadPlanes;

  int _nhit, _trajhit, _truInt;
  int _freeplanes, _occ[1000], _trclusthits[1000];
  double _visEng;
  double _trajpur, _trajEng;
  double _plEng[1000], _trajEngPlan[1000];

  //
  double _detX, _detY, _detZ, _vdetZ, _WLSAtten;
  int _octGeom;

  void set_branches();
  void output_liklihood_info(const vector<cluster*>& hits);
  void traj_like(const vector<cluster*>& hits);
  void out_like();
  //


  ///int _trajs_no;
  double _Xtent; 
  bool _assess_event;///
  ///to get no of trajectories
  // std::vector<Trajectory*> _trajs;///for CA
  std::vector<State> _vPR_seed;

 


};

class chiSorter{
public:
  bool operator()(const Trajectory* T1, const Trajectory* T2){
    if ( T2->quality() > T1->quality() ) return true;
    return false;
  }
};

#endif

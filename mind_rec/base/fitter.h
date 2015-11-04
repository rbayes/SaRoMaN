
/* -*- mode: c++ -*- */
#ifndef _fitter___
#define _fitter___

#include <recpack/RecPackManager.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/LsqFitter.h>
#include <recpack/ParticleState.h>

#include <mind/MINDsetup.h>
#include <mind/Utilities.h>
#include <mind/MINDfitman.h>
#include <mind/cluster.h>
#include <mind/event_classif.h>
#include <mind/hit_clusterer.h>

#include <bhep/event.h>
#include <bhep/gstore.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>

using namespace Recpack;

class fitter{
  
public:
  
  fitter(const bhep::gstore& pstore,
	    bhep::prlevel vlevel=bhep::NORMAL);
    
  virtual ~fitter();
    
  //------------------ main functions ---------------// 
  //bool initialize(const bhep::sstore&) ;
  void Initialize();
  bool Execute(bhep::particle& part,int evNo);
  void Finalize() ;
  //-------------------------------------------------// 
  //Getters.    
  
  std::vector<Trajectory*>&  get_trajs(){return _trajs; }///
  Trajectory& get_hadTrajs() {return _hadTrajs;}
 
  // dict::dictionary <double> getQualityMap(){ if (_reseed_ok) return _traj2->qualitymap(); 
  //else return _traj->qualitymap(); }

  //CA traj
  //  const Trajectory& get_CA_traj(){ 
  //    if (get_classifier().get_int_type()==5) return *_traj; }

  //  std::vector<cluster*> get_CA_meas_vec(){ if (get_classifier().get_int_type()==5) return _meas; }
 
  std::vector<EVector>& get_had_unit(){ return _hadUnit; }
  
  double* get_had_eng(){ return _hadEng; }
  double* get_had_profile() { return _transEdx; }
  std::vector<double>& get_nonMuon_edep() { return _nonMuonEdep; }
  //int get_had_hits() { return _hadmeas.size(); }
  std::vector<int>& get_had_hits() { return _nonMuonHits; }
  int get_fail_event(){ return _failEvent; }///
  bool CheckReseed(){ return _reseed_ok; }
  std::vector<EVector>& get_had_RecDir(){ return _showerDir; }
  double showerVertZ(int i) { return _showerVertZ[i]; }
  double showerNplanes(int i) { return _showerNplanes[i]; }
  double showerXtent(int i) { return _showerXtent[i]; }

  cluster* GetMeas(int num){return _meas[num];}
  std::vector<cluster*> &GetMeasVec(){ return _meas; }
  int GetNMeas(){ return (int)_meas.size();}
  
  MINDsetup GetGeom() { return _geom; }
  //
  cluster* GetMeasurement(bhep::hit& hit);
  
  int GetQ(const Trajectory& traj);

  double GetInitialqP(){ return _initialqP; }

  int* getMuonIndex(){ return _muonindex; }
  //Tempory for likelihoods.
  void set_int_type(const string name){
    get_classifier().set_int_type( name ); }
  //
  //calculate momentum from range.
  //void calculate_len_mom(double len, double *mom);
  
  //recpack manager
  
  RecPackManager& man(){
    return MINDfitman::instance().manager();}

  event_classif& get_classifier(){ return _classify; }
  
  //fit twice

  //void setRefit(bool ok){refit=ok;}
   
protected:
  
  //void resetVirtualPlanes(); 

  //read parameters from store
  void ReadParam();
    
  //seed for fit
  void ComputeSeed(const Trajectory& traj,State& seed, int firsthit=0);
  void ComputeMomFromParabola(const Trajectory& traj, int nplanes, int firsthit, EVector& V);
  void ComputeMomFromRange(const Trajectory& traj, int nplanes, int firsthit, EVector& V);

  //seed error
  void ApplyCovarianceFactor(double factor, EMatrix& C0);
 
  //fit trajectory
  bool FitTrajectory(const State& seed, const int trajno);
  //bool ReseedTrajectory(Trajectory& traj,const int trajno);
  bool ReseedTrajectory(const int trajno);
  //bool fitHadrons();
  //double eng_scale(double visEng);

  // hadron
  void rec_had_energy();
  void rec_had_edep(int j);

  //-------- get traj from event----------//
  bool readTrajectory(const bhep::particle& part);

  // Create a single trajectory with all measurements when the PR is off
  bool CreateSingleTrajectory(Trajectory& traj); 

  // Create measurements from hits
  bool CreateMeasurements(const bhep::particle& part); 


  // Check traj passes cuts for fitting.
  bool CheckValidTraj(const Trajectory& traj);


  ///calcute seed for refit
  void ComputeSeedRefit(const Trajectory& traj, State& seedState);
  //string getPlaneName(bhep::hit);
  //--------------------------------------//
  
  bool CheckQuality(const Trajectory& traj);
    
  void Reset();

protected:

  bhep::gstore _store;
  
  bhep::prlevel _level;
    
  bhep::messenger _m;
    
  MINDsetup _geom;
  
  //counter for virtual planes
  //size_t pnumber;
  
  //Parameters to define fitting method.
  bool _refit; //Do second fit.
  bool _patternRec; //Pattern recognition algorithm required?
  
  int _fitCheck;
  int _min_seed_hits; //Minimum isolated hits required for Prec seed.
  double _min_iso_prop;

  //bit to tell if reseed perfromed
  bool _reseed_ok;
  bool _reseed_called ;
  bool _fitted;

  //------------------ Physics -----------------//
    
  double _X0;
  double _rel_dedx;
  double _rel_dens;
  double _widthI;
  double _widthS;

  int dim; //dimension of model state
  //int meas_dim; //dimension of measurement
  //double _tolerance; //pos. resolution/plane tolerance
  
  State _seedstate;   
  //EVector qoverp;
  double _initialqP;
  
  string _model;  // fit model
  //string kfitter; // kind of fit
  
  //double chi2node_max;
  //int max_outliers;
  double _chi2fit_max;
  double _facRef;

  //Detector name for hit getter. made members 26/11, weird error day.
  string _detect;
  int _highPass;
  int _lowPass;
  double _lowFit1, _lowFit2;
  int _muonindex[2];
  std::vector<double> _nonMuonEdep ;
  std::vector<int> _nonMuonHits ;
  std::vector<EVector> _showerDir;
  std::vector<double> _showerVertZ;
  std::vector<double> _showerNplanes;
  std::vector<double> _showerXtent;

  Trajectory _traj;
  Trajectory _traj2;
  Trajectory _traj3;
  std::vector<Trajectory*> _trajs;  //
  Trajectory _hadTrajs;
  


  std::vector<cluster*> _meas;
  std::vector<cluster*> _hadmeas;

 

  //value set to identify where a traj failed:
  //0=too few hits. 1=too many hits. 2=outside fiducial. 3=no convergence with kink.
  //4=Couldn't find seed for patrec. 5=Failed in pat rec filtering
  int _failType;
  int _failEvent;
  int _intType ;///
  int _pr_count;///
  // int reseed_count;

  //size_t nnodes;

  //Temporary fix (??). vector to store hadron unit dir vec.
  std::vector<EVector> _hadUnit;
  double _hadEng[2];
  double _transEdx[2];
  double _length;
 
  // Stuff relevant for event classification, will be uncommented when needed.
  event_classif _classify;

  // hit clustering.
  hit_clusterer *_clusters;

  bool _doClust;
  //
  std::vector<State> _vPR_seed;
  //-------------- verbosity levels ------------//

  //int vfit,vnav,vmod,vmat,vsim;
  
};

class reverseSorter{
public:
  bool operator()(const cluster* p1, const cluster* p2){
    if (p2->position()[2] < p1->position()[2]) return true;
    return false;
  }

};

class forwardSorter{
public:
  bool operator()(const cluster* p1, const cluster* p2){
    if (p2->position()[2] > p1->position()[2]) return true;
    return false;
  }

};
  
#endif


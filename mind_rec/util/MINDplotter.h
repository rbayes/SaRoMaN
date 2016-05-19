
/* -*- mode: c++ -*- */
#ifndef _mind_plotter___
#define _mind_plotter___

#include <bhep/event.h>
#include <bhep/particle.h>

#include <TFile.h>
#include <TTree.h>

#include <recpack/RecPackManager.h>
//#include <recpack/Measurement.h>
#include <recpack/Ring.h>

#include <mind/fitter.h>
#include <mind/cluster.h>

#include "TVector3.h"



using namespace std;
using namespace Recpack;
using namespace bhep;


/* A Class with functions to plot the results of fits to MIND data */

class MINDplotter{

 public:

  MINDplotter();

  virtual ~MINDplotter();

  //Main functions for best output initialization.
  void Initialize(string outFileName,  bool patRec, bool clust, bhep::prlevel vlevel=bhep::NORMAL);
  void Execute(fitter& Fit, const bhep::event& evt);///
  void Finalize();
//   //
  //To calculate expected vertex position give trajectory and vertex location.
  bool extrap_to_vertex(const Trajectory& traj, 
			const bhep::Point3D& vertexLoc, fitter& fitObj, State& ste, const int trajno);

  //To calculate pos, charge, direction and momentum without extrapolating to vertex //tapasi
  void fill_kinematics(const Trajectory& traj, State& ste, const int trajno);

  /*Requires a vector of the node of interest, corresponding covariance Matrix
    and the event under study to calculate either position or momentum pulls */
  void position_pulls(const int trajno);
  void momentum_pulls(const int trajno);
  void direction_pulls(const int trajno);

  /* Function to find maximum local chi2 in a particular trajectory
     and enter that value in a histogram */
  void max_local_chi2(const Trajectory& traj, const int trajno);

  /* Function to plot stats about pattern recogntion, 1: no clust, 2: clust */
  //   void patternStats1(fitter& Fit,const Trajectory& traj, const int trajno ); ///
  void patternStats2(fitter& Fit);///
  void hitBreakUp(fitter& Fit); ///

  /*Function to record quality of hadron fit*/
  void hadron_direction(fitter& fit);

  /*tru interaction type*/
  void get_tru_int_type(const bhep::event& evt);

  /*get event information*/
  void get_event_info(const bhep::event& evt,fitter& fit);

protected:

  bhep::prlevel _level;
    
  bhep::messenger _m;

  TFile *outFile;

  TTree *statTree;

 
  bool _patR, _clu;
  

private:

  EVector _fittedVert;
  EMatrix _vertMat;
  // EVector _positionStart;
  // TVector3 _positionStart[10];
  //std::vector<vector<double> > _positionStart;
  //std::vector<EVector> _positionStart ;
  //std::vector<EVector> _positionStart ;
  double _positionStart[10][3];
  
  bhep::particle* _truPart;
  int _evNo;
  int _Fit[10];
  int _reFit[10];
  int _fail[10];
  int _ndf[10];
  int _failEvent;
  int _intType[10];
  double _vert[3];
  double _nuEng;
  int _charm[2];
  double _Q2;
  int _pdg[3];
  double _visEng;
  double _engTraj[10];
  double _engvar[4][10];
  //Information on hadrons.
  double _hadP[3];
  double _engTrans;
  int _nhad[2];
  double _chadP[4];
  double _nhadP[4];
  double _hrecP[2][3];
  double _hadE[4][2];
  double _haddot[2];
  double _tX[2][10];
  double _X[2][2][10];
  double _recUnitDir[10][3];
  double _tTh[2][10];
  double _tUnitDir[10][3];
  double _Th[2][2][10];

  double _tqP[10];
  double _qP[3][10];
  
  double _leng[10];
  double _rangqP[2][10];
  int _tQ[10];
  int _tPDG[10];
  int _Q[2][10];

  double _Chi[3][10];
  int _plns[2][10];
  int _nallhits;
  int _nhits[10];
  int _hitType[3][10];
  double _hadEInNucleus;
  

 
  std::vector< vector<double> > _XPos;
  std::vector< vector<double> > _YPos;
  std::vector< vector<double> > _ZPos;
  std::vector< vector<double> > _Edep;
  std::vector< vector<double> > _HTime;
  std::vector< vector<double> > _NodeFitted;

  std::vector<double> _XMeas;
  std::vector<double> _YMeas;
  std::vector<double> _ZMeas;
  std::vector<double> _EMeas;
  std::vector<double> _MuProp;

  std::vector<string> _MotherMeas;
  std::vector<double> _MotherProp;

  double _Theta1;
  double _Theta2;
  double _Theta3;

  double _Delta1;
  double _Delta2;


  std::vector<double> _XHadPos;
  std::vector<double> _YHadPos;
  std::vector<double> _ZHadPos;

 
  int _mus[10][2500];
  bool _cand[10][2500];
  bool _node[10][2500];
  bool _had[10][2500];
  double _pChi[10][3];
  //TString _intName;
  int _truInt;
 
  double _Xtent[10];///
  double _vertZ[10];///
  int _traj_hits;

  
  int _trajs_no; ///
  int _hadHits[10];
  int _hitsInEvent;
  double _nonMuonEdep[10];
  // std::vector<Trajectory*> _trajs;

  int _muontraj[2];

  /*double _chi2[50];
  int _hits[50];
  int _fitted_hits[50];
  int _traj_index;*/
    

  void define_tree_branches();
  //1: Old method. 2: with clustering etc.
  bool extract_true_particle1(const bhep::event& evt, const int trajNo);
  bool extract_true_particle2(const bhep::event& evt);
  void add_to_hads(const bhep::particle& part);
  //
  int _truHitIndex[2][10]; 
  int _hitTrMu;
  int _hitHad;
  double _showerDir[10][3] ;
 
};

#endif


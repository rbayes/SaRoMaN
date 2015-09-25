/*********************************************************
 * Compiled version of the golden channel efficiency and *
 * background analysis of MIND.                          *
 *                                                       *
 * Execute with a parameter file containing file and     *
 * cut information.                                      *
 *                                                       *
 * Andrew Laing, March 2010.                             *
 *                                                       *
 *********************************************************/
/* -*- mode: c++ -*- */

#ifndef _Golden_Analysis__
#define _Golden_Analysis__

#include <bhep/gstore.h>
#include <TRandom3.h>
#include <CLHEP/Random/RanluxEngine.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TMath.h>

#include <vector>

class Golden_Analysis{

 public:
  Golden_Analysis(const bhep::gstore& store);

  ~Golden_Analysis(){};

  void execute();

  void finalize();

 private:
  //functions.
  void define_histos(const bhep::gstore& store);

  void define_rec_binning(const bhep::gstore& store, double* recB, double* truB);

  void define_cut_levels(const bhep::gstore& store);

  void define_variables(const bhep::gstore& store);

  void set_branches(TTree& inTree);

  void analyse_success();

  void fill_tracks2();

  bool fiducial(int ID);

  bool full_cut(int ID);

  bool NC_likelihood(int ID);

  bool kin_cut(int ID);

  bool quad_and_transLong(int ID);
  
  bool include_event();

  //log likelihood values for NC rejection.
  double calc_hitlog();
  void calc_visible_eng();
  void calc_hit_proportion();
  double correctEdep(double edep, double X, double Y);
  double calc_hitfraclog();
  double calc_hitmeanlog();
  double calc_all3log();

  //QT calculation
  double calcQt(double sig);

  //Parabola fit.
  double fit_check(const vector<double>& XP, const vector<double>& ZP);

  //Quasi-elastic formulae.
  double quasi_nubar_prot();
  double quasi_nu_neut();

 private:
  //members.

  //Common base of all files to be read.
  std::string _filebase;

  //Output file.
  TFile* _OutFile;

  //Likelihood files.
  TFile* _likeCC;
  TFile* _likeNC;

  //Analysis Type: 1 - hit only in NC like, 2 - all NC like.
  int _anType;

  //Out Histograms.
  TH1F *_Eff;
  TH1F *_Ey;
  TH1F *_rcuts;
  TH1F* _pcuts;
  TH1F *_rcuts_sig;
  TH1F* _pcuts_sig;
  TH1F *_rcuts_back;
  TH1F* _pcuts_back;
  TH1F* _ifail;
  TH1F* _fitprop;
  TH1F* _zFid;
  TH1F* _zFidcut;
  TH1F* _Lqp;
  TH1F* _Lqp_back;
  TH1F* _LCC1;
  TH1F* _LCC2;
  TH1F* _LCC3;
  TH1F* _LCC4_sig;
  TH1F* _LCC4_back;
  TH1F *_back;
  TH2F *_back_pion;
  TH2F *_back_kaon; 
  TH1F *_By;
  TH2F *_recEff;
  TH2F *_recBack;
  TH1F *_pmu;
  TH1F *_pmufail;
  TH1F *_pmutfail;
  TH2F *_recmunu;
  TH2F *_recQtnu;

  TH2D* _dir_back;
  TH2D* _dir_sig;
  TH2D* _dirMC_back;
  TH2D* _dirMC_sig;

  TH1F *_alltracks1;
  TH2F *_alltracks2;
  TH1F *_ally;
  TH1F *_qCharge;
  TH1F *_qCharge_back;
  //
  //Likelihood histograms.
  TH1F *_hitCC, *_hitNC;
  TH1F *_fracCC, *_fracNC;
  TH1F *_meanCC, *_meanNC;
  TH1F *_errSig, *_errBack;
  //
  std::vector<TH1F*> IntHist;
  std::vector<TH1F*> IntEff;
  std::vector<TH1F*> Intback;

  TH1F* _EffNM;
  TH1F* _alltracksNM1;
  TH1F* _backNM;

  //Interaction type: nu_mu, NC, nu_e.
  std::string _intType;

  //Charge expected for beam nu_mu interaction.
  int _beamCharge;

  //Had rec or smear: 0, 1 (resp.).
  bool _doSmear;
  
  //Random Generator for the smear.
  TRandom3 _ranGen;

  //Max energy of neutrinos.
  double _maxE;
  double _maxOver;

  //First and last file number.
  int _firstF;
  int _lastF;

  //Detector parameters.
  double _detX;
  double _detY;
  double _WLSAtten;

  //Stuff for the xsec systematic study.
  int _doVary; //If you want it to be done(0), reduce Vtype (1), reduce 'rest' (2).
  int _Vtype; //Which interaction is of interest.
  TFile *_xsecF;
  TH1F *_xsecErr;//The binned errors.
  //

  //Variables to be read from the tree.
  char _suc, _inCand[2500], _hadH[2500];
  int _charge[3], _trHits[5], _nhits, _truInt, _fail;
  double _enu, _mom[3], _hadEng[2], _hadMom[3], _maxLoc[2], _dir[3][2];
  Double_t _xPos[2500], _yPos[2500], _zPos[2500], _depEng[2500];
  int _npi[3], _nk[3];
  //
  double _smearRoot, _smearConst;
  double _angA, _angB;//Smear constants for angular smear.
  //Cut levels.
  vector<double> _edges;//for edge/fiducial cut.
  double _logQPmin;//Minimum err/mom log likelihood value .
  double _minNodes;//minimum fraction of nodes fitted.
  vector<double> _likeCuts;//Min NC logLike values.
  double _kMinE;//min eng for kin cuts.
  double _minQT;//minimum Qt variable.
  double _momGrad;//Gradient for calculation of minmom.
  vector<double> _quadFitCuts;//must be <[0] or >[1].
  vector<double> _translongCuts;//if hits<[0], trans/long>[1].
  //
  //Other parameters to be calculated.
  double _trajEng;
  double _visEng;
  double _trHitProp;
  double _recEnu;

  bool isMeson;
};

#endif

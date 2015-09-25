
#ifndef _GoldCuts_TrajSel__
#define _GoldCuts_TrajSel__

#include <bhep/gstore.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TMath.h>

#include <vector>
#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;

class GoldCuts_TrajSel{

 public:
  GoldCuts_TrajSel(const bhep::gstore& store);
  ~GoldCuts_TrajSel(){};
  void execute();
  void finalize();

 private:
  
  void define_histos(const bhep::gstore& store);
  void define_rec_binning(const bhep::gstore& store, double* recB, double* truB);
  bool include_event();
  void define_cut_levels(const bhep::gstore& store);
  void define_variables(const bhep::gstore& store);
  void set_branches(TTree& inTree);
  void select_muon_trajectory();
  void analyse_success();
  void fill_tracks2();
  bool fiducial(bool usetrue);
  bool full_cut();
  // bool Muon_Selection();
  bool NC_likelihood(int ID, double& Pratio);
  double calc_hitlog();
  bool kin_cut(int ID);
  double quasi_nubar_prot();
  double quasi_nu_neut();
  double calcQt(double sig);
  //members.

 private:
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
  TH2F* _imuontraj;
  TH2F* _itraj;
  TH1F* _fitprop;
  TH1F* _fitprop_back;
  TH1F* _zFid;
  TH1F* _zFidcut;
  TH2F* _rzFid;
  TH2F* _rzFidcut;
  TH1F* _Lqp;
  TH1F* _Lqp_back;
  TH1F* _errqp;
  TH1F* _errqp_sig;
  TH1F* _errqp_back;
  TH1F* _LCC1_sig;
  TH1F* _LCC1_back;
  TH1F* _Nhits;
  TH1F* _Nhits_sig;
  TH1F* _Nhits_back;
  TH2F* _rangtruE;
  /*TH1F* _LCC2_sig;
  TH1F* _LCC2_back;
  TH1F* _LCC3_sig;
  TH1F* _LCC3_back;
  TH1F* _LCC4_sig;
  TH1F* _LCC4_back;*/
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
  TH2F *_recmutrunu;
  TH2F *_recQtnu;
  TH2F *_recmunu_back;
  TH2F *_recQtnu_back;

  TH2D* _dir_back;
  TH2D* _dir_sig;
  TH2D* _dirMC_back;
  TH2D* _dirMC_sig;

  TH1F *_alltracks1;
  TH2F *_alltracks2;
  TH1F *_ally;
  /*TH2F* _dispX_sig;
  TH2F* _dispZ_sig;
  TH2F* _dispX_back;
  TH2F* _dispZ_back;
  TH1F *_qCharge;
  TH1F *_qCharge_back;*/
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
  std::vector<TH2F*> IntrecEff;
  std::vector<TH2F*> IntrecBack;

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
  // char _; //  _inCand[2500], _hadH[2500];
  int _evt;
  int _suc[10];
  double _enu;
  int _ntraj;
  int _mutrajlong;
  int _mutrajNhit;
  int _mutraj;
  int _truq;
  int _recq[10];
  int _truID[10];
  int _truInt;
  int _fail;
  double _truqP;
  double _recqP[10];
  double _errqP[10];
  double _rangeMom[10];
  double _initqP[10];
  double _recXdir[10];
  double _recYdir[10];
  double _truXdir;
  double _truYdir;
  int _fitn[10];
  int _inmu[10];
  int _hadhits;
  double _hadTruEng;
  double _hadRecEng[2];
  double _hadQESEng[2];
  double _hadEng;
  double _hadMom[3];
  double _hrdir[3];
  int    _nohits[10];
  std::vector<vector<double> >* _XPos;
  std::vector<vector<double> >* _YPos;
  std::vector<vector<double> >* _ZPos;
  double _len[10];
  double _vertX[10];
  double _vertY[10];
  double _vertZ[10];
  double _trueVert[3];
  //
  double _smearRoot, _smearConst;
  double _angA, _angB;//Smear constants for angular smear.
  //Cut levels.
  vector<double> _edges;//for edge/fiducial cut.
  double _logQPmin;//Minimum err/mom log likelihood value .
  double _minNodes;//minimum fraction of nodes fitted.
  double _rangelimit;
  double _likeCut;//Min NC logLike values.
  double _kMinE;//min eng for kin cuts.
  double _minQT;//minimum Qt variable.
  double _momGrad;//Gradient for calculation of minmom.
  vector<double> _quadFitCuts;//must be <[0] or >[1].
  vector<double> _translongCuts;//if hits<[0], trans/long>[1].
  //
  // Calculated variables.
  double QT, _recEnu;
  std::vector<int> traj;
  std::vector<int>::iterator itraj;

};

#endif

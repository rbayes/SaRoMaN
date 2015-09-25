/**************************************************************
 *                                                            *
 * Class to perform smear and cut based analysis on MIND      *
 * data.                                                      *
 *                                                            *
 * Authors: A. Laing, A. Cervera; February 2009.              *
 *                                                            *
 **************************************************************/

/* -*- mode: c++ -*- */
#ifndef _old_analysis___
#define _old_analysis___

#include <bhep/event.h>
#include <bhep/particle.h>
#include <bhep/gstore.h>
#include <bhep/messenger.h>

#include <TRandom3.h>

#include <TFile.h>
#include <TTree.h>

#define Pi 3.1415926535897932

class old_analysis{

 public:

  old_analysis( bhep::prlevel level);

  virtual ~old_analysis();

  //--------------- main functions -----------------//
  bool initialize(const bhep::gstore& rstore,
		  const bhep::gstore& astore);
  bool execute(const bhep::event& evt);
  bool finalize();
  //

 protected:

  TRandom3 ranGen;

  bhep::gstore smear_store;

  bhep::gstore cut_store;

  bhep::messenger m;

  //Output file/tree
  TFile* outFile;
  TTree* dataOut;

  void readParam();
  void define_tree_branches();

  void get_event_properties(const bhep::event& evt);
  void get_particle_info();
  double check_daughters(const vector<const bhep::particle*>& hijas);
  void rec_Energy(const bhep::event& evt);
  double rec_had_eng(double eng);
  double rec_mu_eng();
  void smear_ang(bhep::Vector3D& muVec, bhep::Vector3D& hadVec);
  void rec_mu_angle(bhep::Vector3D& muVec);
  void rec_had_angle(bhep::Vector3D& hadVec);
  int get_int_type(string intName);
  bool cut_length();
  bool kin_cuts();

  //members from parameters.
  double _hadSmearSig1;
  double _hadSmearSig2;
  double _X0;
  double _measRes;
  double _lambFe;
  double _had_ang_sig[2];
  double _Bfield;

  //cuts.
  double _minlenDiff;
  double _minQt;
  double _momGrad;
  //double _sigCharge;

  //Event parameters used more than once.
  double _leadL;
  double _leadLproj;
  int _nleadhits;
  double _cosAng;
  double _HadEng;
  double _muhadang;
  double _truMom;
  double _muEng;
  int _dHits;
  vector<bhep::particle*> _parts;
  vector<bhep::particle*>::iterator _LpIt;

  //output parameters.
  int _intType;
  bool _CC;
  double _recEng[2];
  double _leadp[2];
  double _Qt;
  double _lendiff;

};

#endif

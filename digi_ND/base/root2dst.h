/* -*- mode: c++ -*- */
#ifndef _root_dst___
#define _root_dst___

#include <bhep/event.h>
#include <bhep/particle.h>
#include <bhep/messenger.h>
#include <bhep/bprint.h>
#include <bhep/bhep_svc.h>
#include <bhep/ray.h>
#include <bhep/hit.h>
#include <bhep/gstore.h>
//#include <bhep/clhep.h>
//#include <CLHEP/Random/RanluxEngine.h>

#include <digi/gdml_hit_constructor.h>
//#include <digi/hit_constructor.h>

#include <TTree.h>
#include <Riostream.h>
#include <sstream>

#include "TFile.h"
#include "TH1F.h"

using namespace bhep;

//! root2dst Class
/*!
  Converts ROOT tree into bhep DST file
*/

class root2dst{
  
public:
  
  root2dst(bhep::prlevel, bhep::gstore* store=NULL);
  
  ~root2dst(){};
  
  bool initialize(double *res, long seed, TTree *InPutTree=NULL,
		  TString OutFileName="");
  bool execute();
  bool finalize();

  //Functions to create bHEP object from the ntuple
  void createEvent();
  bool make_particles();
  particle* define_lead_particle();
  bool hits_fromFile(vector<hit*>& muHit, vector<hit*>& hadHit);
  particle* define_hadron();
  particle* create_digital_representation(const vector<particle*>& tru_parts);

  void print();

  TFile* digiOutfile;
  //TList* hList;
  TH1F* rawHits;
  TH1F* clusteredHits;
  TH1F* digitizedHits;
  TH1F* xeTH1F;
  TH1F* xeAttTH1F;
  TH1F* xeSmearTH1F;
  TH1F* yeTH1F;
  TH1F* yeAttTH1F;
  TH1F* yeSmearTH1F;
  vector<TH1F*> histo_vec;



  //Handy int/float to string converters.
  // TString ToString(Int_t num){
//     ostringstream start;
//     start<<num;
//     TString start1=start.str();
//     return start1;
    
//   }
  
//  TString ToString(Float_t num){
//     ostringstream start;
//     start<<num;
//     TString start1=start.str();
//     return start1;
    
//   }
  
protected:
  
  //verbosity level
  bhep::prlevel level;
  
  //messenger
  bhep::messenger m;
  
  //event counter
  size_t nevt;
  
private:

  //The Event
  bhep::event nuEvent;

  //file to write dst to.
  //writer_gz outgz;

  //Root tree to be read.
  TTree *dataIn;

  //Random engine for smearing.
  //RanluxEngine ranGen;
  //Voxel hit constructor for digitization.
  gdml_hit_constructor* _construct;
  //hit_constructor* _construct;
  //Are voxels to be made.
  bool _doVox;

  //Gaussian sigma for smearing position.
  double sigMa;
  //Gaussian sigma for smearing deposited energy.
  double sigMaE;
  //hit check
  Int_t NmuHits, NhadHits;

};





#endif

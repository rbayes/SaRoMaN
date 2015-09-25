
#ifndef _Train_TMVA__
#define _Train_TMVA__

#include <bhep/gstore.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TObjString.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCut.h>
#include <TStopwatch.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodCuts.h"

#include <vector>

class Train_TMVA{

 public:
  Train_TMVA(const bhep::gstore& store);
  ~Train_TMVA() {};
  void execute();
  void finalize();

 private:

  void set_MVA_methods(const bhep::gstore& store);
  void define_cuts(const bhep::gstore& store);
  
 private:
  TFile* _OutFile;
  std::string _filebaseSig;
  int _firstSig;
  int _lastSig;
  std::vector<std::string> _filebaseBack;
  std::vector<int> _firstBack;
  std::vector<int> _lastBack;
  std::map<std::string,int> Use;

  // Data manipulation
  TMVA::Factory* _factory;
  
  // Parameters for cuts
  std::string _intType;
  int         _anType;
  int         _beamCharge;
  double      _maxE;
  double      _maxOver;
  double      _minNodes;
  
  // A set of fixed cuts on the analysis.
  TCut _cuts;
};

#endif

#ifndef _QGSP_BERT__
#define _QGSP_BERT__

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

class MindQGSP_BERT: public G4VModularPhysicsList
{

 public:
  MindQGSP_BERT(G4int ver = 1);
  virtual ~MindQGSP_BERT();

 public:

  void SetCuts();
  void ConstructProcess();

 private:

  void ConstructUserLimits();

};

#endif

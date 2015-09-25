/* -*- mode: c++ -*- */
#ifndef _setupsk___
#define _setupsk___

#include <recpack/RecPackManager.h>
#include <bhep/messenger.h>
#include <bhep/gstore.h>

using namespace Recpack;

class NSetupSk{

public:

  virtual ~NSetupSk(){};
    
  virtual Setup& setup() = 0;
    
protected:

  Setup _gsetup;
    
  bhep::messenger _msetup;
    
  //store for geom params
  //bhep::sstore _store;

};



#endif 

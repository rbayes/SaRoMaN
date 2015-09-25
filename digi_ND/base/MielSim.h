/* -*- mode: c++ -*- */
#ifndef _MielSim___
#define _MielSim___

#include <bhep/event.h>
#include <bhep/particle.h>
#include <bhep/messenger.h>
#include <bhep/bprint.h>
#include <bhep/bhep_svc.h>
#include <bhep/ray.h>
#include <bhep/hit.h>
//#include <bhep/clhep.h>
#include <TRandom3.h>
#include <bhep/gstore.h>

//#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <sstream>

using namespace bhep;

//! MielSim Class
/*!
  bla bla blaaaa....
*/

class MielSim {
  
public:
  
  MielSim(bhep::prlevel);
  ~MielSim(){};
  
  bool Initialize(const bhep::gstore& run_props);
  bhep::particle* CreateDigitalRepresentation(const vector<bhep::particle*>& tru_parts);
  bool Finalize();

protected:
  //verbosity level
  bhep::prlevel fVerbosityLevel;
  //messenger
  bhep::messenger fMessenger;
  //event counter
  size_t fEventNumber;
  //simulation of the electronic responce 
  std::vector<bhep::hit*> DigitizeHits(const std::vector<bhep::hit*>& trueHits);
  
private:
  //The Events
  vector<bhep::event*> nuEvent;
  //File to write bhep-dst to
  //bhep::writer_gz fZippedOutPutFile;
  //Random engine for smearing.
  TRandom3 fRandomGenerator;
  // The scintillator bar simulation 
  //ScintBar* fScintBar;
//   // The sensor simulation 
//   PMSensor* fPMSensor;
//   // The electronics simulation
//   Electronics* fElectronics;
};





#endif

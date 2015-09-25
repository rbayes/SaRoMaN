/* -*- mode: c++ -*- */
#ifndef _ScintBar___
#define _ScintBar___

#include <bhep/hit.h>
#include <bhep/clhep.h>
#include <bhep/messenger.h>
#include <TString.h>
#include <Riostream.h>
#include <sstream>

//using namespace bhep;

//! ScintBar Class
/*!
  bla bla blaaaa....
*/

class ScintBar {     //Should we here inherit from a base class to make it possible 
                     //to have same structure for INO?
       
public:
  
  ScintBar(bhep::prlevel verbosityLevel);
  ~ScintBar(){};
  
  bool Initialize();
  std::vector<WLSPhoton*> HitsToWLSPhotons(std::vector<bhep::hit>& hits);

private:
  // Verbosity level
  bhep::prlevel fVerbosityLevel;
  // Messenger
  bhep::messenger fMessenger;
  // Average number of photons trapped in WLS bar for 1 MeV/cm energy deposit 
  double fPhotPerMeV;
  // Reflection Probability
  double fReflectionProb;
  // The effective light velocity in the scintillator and WLS combo.
  double fLightVelocity;
  // The Birks' law coefficient.
  double fBirksConstant;
  // The WLS fibre decay time
  double fDecayTime;
  // The Scintilator Bar Dimension
  TVector3 fBarDim;  
  // The long attenuation component of the exponential attenuation model.
  double fLongComponent;   
  // The fraction of light associated with long attenuation component.
  double fLongFraction;
  // The short attenuation component of the exponential attenuation model.
  double fShortComponent;
  // The spline to save the attenuation curve.
  TSpline3* fAttenuation;

};





#endif

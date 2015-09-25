#include <ScintBar.h>

using namespace std;
using namespace bhep;

//*************************************************************
ScintBar::ScintBar(bhep::prlevel verbosityLevel) {
  fVerbosityLevel = verbosityLevel;
  fMessenger = bhep::messenger(fVerbosityLevel);
  fMessenger.message("+++ ScintBar Constructor +++", bhep::VERBOSE);
}

//*************************************************************
bool ScintBar::Initialize() { 
  fMessenger.message("+++ ScintBar init  function ++++", bhep::NORMAL);
  //Should read in from some DB
  fPhotPerMeV = 25;
  fReflectionProb = 0.99;
  fLightVelocity = 20*bhep::cm/bhep::ns;
  fBirksConstant = 0.005*bhep::cm/bhep::MeV;
  fDecayTime = 11*bhep::ns;
  fBarDim.SetXYZ(5*bhep::mm,5*unit::mm,50*unit::cm);  
  fLongComponent = 463.4*bhep::cm;   
  fLongFraction = 0.77;
  fShortComponent = 33.2*bhep::cm;
  fAttenuation = NULL;
  return true;
}

//*************************************************************
vector<WLSPhoton*> HitsToWLSPhotons(vector<bhep::hit>& hits) {
  //Loop over the hits
  for (vector<bhep::hit>::iterator hitIt = hits.begin(); hitIt != hits.end(); hitIt++) {
    const bhep::Point3D hitPoint = hitIt->x();
    
  }




}

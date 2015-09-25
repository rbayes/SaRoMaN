#include <MielSim.h>
#include <TRandom3.h>

using namespace std;
using namespace bhep;

//*************************************************************
MielSim::MielSim(bhep::prlevel verbosityLevel) {

  fVerbosityLevel = verbosityLevel;
  fMessenger = bhep::messenger(fVerbosityLevel);
  fMessenger.message("+++ MielSim Constructor +++", bhep::VERBOSE);
  
  //ScintBar* fScintBar = new ScintBar(bhep::prlevel verbosityLevel);
  //PMSensor* fPMSensor = new PMSensor(bhep::prlevel verbosityLevel);
  //Electronics* fElectronics = new Electronics(bhep::prlevel verbosityLevel);
}

//*************************************************************
bool MielSim::Initialize(const bhep::gstore& run_props) { 
  fMessenger.message("+++ MielSim init  function ++++", bhep::NORMAL);
  //This function will initiate all the relevant parameters via the bhep store
  long seed;
  if ( run_props.find_dstore("Gen_seed") )
    seed = (long)run_props.fetch_dstore("Gen_seed");
  fRandomGenerator = TRandom3(seed);
  //Initialize event number
  fEventNumber=0;
  return true;
}

//*************************************************************
particle* MielSim::CreateDigitalRepresentation(const vector<particle*>& tru_parts) {
  //Create the new particle that is going to be the digital 
  //representation of the input particle
  ptype pT = DIGI;
  particle *digiParticle = new particle(pT, "unknown");

  vector<particle*>::const_iterator partIt;
  vector<bhep::hit*>::const_iterator digiHitIt;

  for (partIt = tru_parts.begin(); partIt != tru_parts.end();partIt++) {
    //Get the true hits of the particle
    //const vector<bhep::hit*>& trueHits = (*partIt)->hits("MIND");
    //Digitize, i.e. get the digitized hits of the particle
    vector<bhep::hit*> digiHits = DigitizeHits( (*partIt)->hits("MIND") );
    //Put the digitized hits in the digitized particle
    for (digiHitIt = digiHits.begin(); digiHitIt != digiHits.end(); digiHitIt++) {
      (*digiHitIt)->set_mother_particle( *(*partIt) );
      digiParticle->add_hit("MIND", (*digiHitIt));
    }
  }
  
  return digiParticle;
}

//*************************************************************
vector<bhep::hit*> MielSim::DigitizeHits(const vector<bhep::hit*>& trueHits) {

  vector<bhep::hit*> digiHits;
  //TEMP
  //Loop over the true hits 
  vector<bhep::hit*>::const_iterator hitIt;
  for (hitIt = trueHits.begin(); hitIt != trueHits.end(); hitIt++) {
    double sigma = 1.0;  //TEMP anyway
    double X = (*hitIt)->x().x()/cm + fRandomGenerator.Gaus(0, sigma);
    double Y = (*hitIt)->x().y()/cm + fRandomGenerator.Gaus(0, sigma);
    double Z = (*hitIt)->x().z()/cm;
    Point3D hitPos(X * cm, Y * cm, Z * cm);
    hit* digitizedHit = new hit("MIND");
    digitizedHit->set_point(hitPos);
    digiHits.push_back(digitizedHit);
  }
  return digiHits;  //Don't know any reason to fail... yet...
  //END TEMP

  // Convert G4 hits (bhep::hit) to WLS photons
  //vector<WLSPhoton*> wlsPhotons = fScintBar->HitsToWLSPhotons(trueHits);
   
  // Convert WLS photons to avalanches in the photo multiplier sensor
  //vector<QTSignals*> sensorQTSignals = fPMSensor->GenerateAvalanches(wlsPhotons);
 
  // Convert the photo multiplier charge time signals (i.e. avalanches) to 
  // electronic charge time signals (i.e. electronic pulses)
  //vector<QTSignals*> elecQTSignals = fElectronics->GeneratePulses(sensorQTSignals);

  /// ????? 
  //Should MielSim produce only electronic signals or should there be one more
  //step here where the electronic signals are de-digitized, i.e. should MielSim
  //reconstruct the digitized hits? 

}

//***************************************************************************
bool MielSim::Finalize() {

  return true;
}

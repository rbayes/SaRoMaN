/* Modified from QGSP_BERT.icc
   Andrew Laing, 2009 */

#include "MindQGSP_BERT.h"
#include "MindConfigService.h"
#include "MindParamStore.h"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4DataQuestionaire.hh"
// #include "HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4UserSpecialCuts.hh"
//
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

MindQGSP_BERT::MindQGSP_BERT(G4int ver) : G4VModularPhysicsList()
{
  
  G4DataQuestionaire it(photon);
  G4cout << "<<< Geant4 Physics List simulation engine: QGSP_BERT 3.3"<<G4endl;
  G4cout <<G4endl;
  
  this->defaultCutValue = 0.7*mm;  
  this->SetVerboseLevel(ver);
  
  // EM Physics
  this->RegisterPhysics( new G4EmStandardPhysics(ver));
  
  // Synchroton Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics("extra EM"));
  
  // Decays
  this->RegisterPhysics( new G4DecayPhysics("decay",ver) );
  
  // Hadron Elastic scattering
  this-> RegisterPhysics( new G4HadronElasticPhysics(ver));
  
  // Hadron Physics
  G4bool quasiElastic;
  // this->RegisterPhysics( new HadronPhysicsQGSP_BERT("hadron",quasiElastic=true));
  this->RegisterPhysics( new G4HadronPhysicsFTFP_BERT("hadron",quasiElastic=true));

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics("stopping"));
  
  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics("ion"));
  
  // Neutron tracking cut
  this->RegisterPhysics( new G4NeutronTrackingCut("Neutron tracking cut", ver));
  
}

MindQGSP_BERT::~MindQGSP_BERT()
{
}

void MindQGSP_BERT::SetCuts()
{
  if (this->verboseLevel >1){
    G4cout << "QGSP_BERT::SetCuts:";
  }
  
  //this->SetCutsWithDefault();   
  const MindParamStore& config = MindConfigService::Instance().Physics();
  G4double prodCut = config.GetDParam("production_cut") * mm;

  SetCutValue(prodCut, "gamma");
  SetCutValue(prodCut, "e-");
  SetCutValue(prodCut, "e+");
  
  // if (this->verboseLevel >0)
//     G4VUserPhysicsList::DumpCutValuesTable();  
}

void MindQGSP_BERT::ConstructProcess()
{
  AddTransportation();
  /*
  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    (*itr)->ConstructProcess();
  }
  */
  ConstructUserLimits();

}

void MindQGSP_BERT::ConstructUserLimits()
{

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (particle->GetParticleName()!="gamma" &&
	particle->GetParticleName()!="mu+" &&
	particle->GetParticleName()!="mu-")
      pmanager->AddProcess(new G4UserSpecialCuts,     -1,-1,1);

  }
}

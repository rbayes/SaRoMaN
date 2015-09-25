#ifndef __MIND_CRY_INTERFACE__
#define __MIND_CRY_INTERFACE__ 1

#include <G4VPrimaryGenerator.hh>
#include <globals.hh>

#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
// #include "Utils/CmdLineArgParserUtils.h"
// #include "Utils/CmdLineArgParserException.h"
#include "Utils/CmdLnArgParser.h"

#include <bhep/bhep_svc.h>

class G4ParticleDefinition;
class G4PrimaryParticle;
class G4Event;

class MindCryInterface: public G4VPrimaryGenerator
{
 public:
  // Constructor
  MindCryInterface();
  // Destructor
  ~MindCryInterface();

  void GeneratePrimaryVertex(G4Event*);
  
 private:
  void Initialize();

  G4int SelectVertexRegion();
  G4double GetTargetProb();

  G4ThreeVector GenerateRandomDirection();
  G4double GenerateRandomEnergy(G4double, G4double);

 private:
  std::vector<CRYParticle*> *vect;
  G4ParticleTable* particle_table;
  CRYGenerator::CRYGenerator* gen;

  G4ParticleDefinition* _particle_definition;
  G4double _mass; ///< Particle mass
  G4double _charge; ///< Particle charge
  G4int _vert_mat;
  G4double _FeTargetProp;
  G4String _vtx_location;

  G4int _evCount[2]; //< Event counter.
  G4int _surface;
  // G4int _tasd; //< a tasd detector?

  //Vertex in the case a fixed point is requested.
  G4int _fvecReg;
  G4ThreeVector _fvec;
  bhep::vdouble _fspdg;
};

#endif

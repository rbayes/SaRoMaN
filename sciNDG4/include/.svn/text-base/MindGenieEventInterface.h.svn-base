#ifndef __GENIE_EVENT_INTERFACE__
#define __GENIE_EVENT_INTERFACE__

#include <G4VPrimaryGenerator.hh>
#include <globals.hh>
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4Material.hh>
#include <G4PrimaryParticle.hh>

#include "TGeoManager.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"

#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "Utils/XSecSplineList.h"
#include "Utils/UnitUtils.h"
#include "Utils/StringUtils.h"
#include "Geo/PointGeomAnalyzer.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "EVGDrivers/GeomAnalyzerI.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "EVGCore/EventRecord.h"
#include "FluxDrivers/GSimpleNtpFlux.h"
#include "Messenger/Messenger.h"

#include "MindConfigService.h"
#include "MindParamStore.h"

#include <bhep/bhep_svc.h>

class G4ParticleDefinition;
class G4PrimaryParticle;
class G4Event;

using namespace genie;
// using namespace genie::controls;

class MindGenieEventInterface : public G4VPrimaryGenerator 
{
 public:
  MindGenieEventInterface();
  virtual ~MindGenieEventInterface();

  void GeneratePrimaryVertex(G4Event* event);

 private:
  void Initialize();
  void Finalize();

  void initRandGen();
  void initGenieMessenger();
  void initGenieSplines();
  void initGenieGeomDriver();
  void initGenieFluxDriver();
  void initGenieMCjobDriver();
  void initGenieMCjobMonitor();

  G4PrimaryParticle*
    CreatePrimaryParticle(GHepParticle& part, G4int PDG, bool primLep=false);
  
 private:
  std::string _geom_root_file;
  std::string _musr_flux_tree;
  std::string _genie_splines;
  bhep::vdouble _had4P;
  bhep::vdouble _fspdg;

  RandomGen*                  _randGen;
  XSecSplineList*             _xspl;
  GeomAnalyzerI*              _geom_anal;
  geometry::ROOTGeomAnalyzer* _root_geom;
  TGeoManager*                _geoMgr;
  // genie::PDGCodeList          _fluxParticles             
  
  flux::GSimpleNtpFlux*       _ntuple_flux_driver;
  
  GFluxI*                     _flux_driver;
  GMCJDriver*                 _mcj_driver;
  GMCJMonitor*                _mcjmonitor;
  EventRecord*                _genie_event;

  double _max_energy;

  G4int _ievent;
  G4int _interactionCount;
  G4int _runNumber;
  G4int _verbose;
  G4int _randInt;

  G4ParticleDefinition* _particle_definition;
  G4double _mass;
  G4double _charge;
  
  G4double _hadEInNucleus;
};

#endif

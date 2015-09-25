
#include "MindGenieEventInterface.h"
#include "MindLookupTable.h"


MindGenieEventInterface::MindGenieEventInterface(){
  Initialize();
}
MindGenieEventInterface::~MindGenieEventInterface(){
  Finalize();
}
void MindGenieEventInterface::Initialize(){
  const MindParamStore& config = MindConfigService::Instance().Generation();
  _runNumber      = config.GetIParam("runNumber");
  _verbose        = config.GetIParam("eventVerbose");
  _max_energy     = config.GetDParam("maxEnergy");

  _geom_root_file = config.GetSParam("geomRootFile");
  _musr_flux_tree = config.GetSParam("musrFluxTree");
  _genie_splines  = config.GetSParam("genieSplines");

  initRandGen();
  initGenieMessenger();
  initGenieSplines();
  initGenieGeomDriver();
  initGenieFluxDriver();
  
  _had4P = bhep::vdouble(4);
  _fspdg = bhep::vdouble(6);
  for(int i=0; i<6; i++) _fspdg[i] = 0;
}
void MindGenieEventInterface::Finalize(){
  // if ( _randGen )           delete _randGen;
  if ( _geom_anal )         delete _geom_anal;
  if ( _root_geom )         delete _root_geom;
  if ( _geoMgr )            delete _geoMgr;
  if ( _ntuple_flux_driver) delete _ntuple_flux_driver;
  if ( _flux_driver )       delete _flux_driver;
  if ( _mcj_driver )        delete _mcj_driver;
  if ( _mcjmonitor )        delete _mcjmonitor;
  if ( _genie_event )       delete _genie_event;
}
void MindGenieEventInterface::GeneratePrimaryVertex(G4Event* event){
  //
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();
  // Get the event
  _genie_event = _mcj_driver->GenerateEvent();
  // Event information
  Interaction* gInt = _genie_event->Summary();
  GHepParticle* fsl = _genie_event->FinalStatePrimaryLepton();
  bevt.add_property( "IntType", gInt->ProcInfo().AsString() );
  
  // In absence of a vertex item in the genie event record, the vertex from the 
  // final state primary lepton is used.
  G4ThreeVector position = G4ThreeVector(fsl->Vx(),fsl->Vy(),fsl->Vz());
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, fsl->Vt());
  // Set Vertex to bhep event.
  bevt.set_vertex( fsl->Vx(), fsl->Vy(), fsl->Vz() );
  // Legacy information from old MIND studies
  bool charmProd = gInt->ExclTag().IsCharmEvent();
  bevt.add_property( "Charm", charmProd );
  if ( charmProd )
    bevt.add_property( "CharmHad", gInt->ExclTag().CharmHadronPdg() );
  // Get the transverse energy exchange Q2
  G4double Q2 = gInt->Kine().Q2(true);
  if ( Q2 != -99999 ) // Protects agains cases where Q2 is not set
    bevt.add_property( "Q2", gInt->Kine().Q2(true)*GeV );
  bevt.add_property( "EngTrans", (_genie_event->Probe()->E() - fsl->E())*GeV );
  //
  
  //Loop over particles.
  TObjArrayIter iter(_genie_event);
  GHepParticle *part = dynamic_cast<GHepParticle *>(iter.Next());//0;
  int pStatus, pdg, fslcount=0;
  _hadEInNucleus = 0;

  while ( part != NULL ){
    //Get particle status.
    pStatus = part->Status();
    
    pdg = part->Pdg();

    if ( pStatus == 0 ) { //initial state particle.
      
      if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) {
	bevt.add_property("nuType",   pdg);
	bevt.add_property("nuEnergy", part->E()*GeV);
      }
      else {
	//Electron scattering events have neutrino, electron and nucleus here.
	if ( pdg >= 1000000000 ) {//Nucleus
	  bevt.add_property("nucType",   pdg);
	  bevt.add_property("nucEnergy", part->E()*GeV);
	} else {
	  bevt.add_property("intpart",   pdg);
	  bevt.add_property("partEnergy", part->E()*GeV);
	}

      }
      
    } else if ( pStatus == 1 && pdg < 2000000000 ) { //stable final state.
      if(fslcount<6)_fspdg[fslcount] = pdg;
      fslcount++;
      if ( part == fsl )
	vertex->SetPrimary(CreatePrimaryParticle( *part, pdg, true ));
      else {
	vertex->SetPrimary(CreatePrimaryParticle( *part, pdg ));
      }
    }
   
    ///hadron energy inside nucleus 
    if(pStatus==kIStHadronInTheNucleus) 
      _hadEInNucleus += part->E();
    
    
    part = dynamic_cast<GHepParticle *>(iter.Next());
  }
  bevt.add_property("fspdg", _fspdg);
  bevt.add_property("had4vec", _had4P);
  bevt.add_property("hadEInNucleus", _hadEInNucleus);
  
  event->AddPrimaryVertex(vertex);
  
  _interactionCount++;
  _genie_event->Clear();

}
void MindGenieEventInterface::initRandGen(){
  _randGen = RandomGen::Instance(); // get random generator initialized
  _randInt = rand() % 100000 + 1;   // generate random number for seed
}

void MindGenieEventInterface::initGenieSplines(){
  _xspl = XSecSplineList::Instance();
  _xspl->LoadFromXml(_genie_splines);

}

void MindGenieEventInterface::initGenieGeomDriver(){
  
  _geoMgr = new TGeoManager();
  _geoMgr->Import(_geom_root_file.c_str());
  _root_geom = new geometry::ROOTGeomAnalyzer(_geoMgr);
  _geom_anal = dynamic_cast<GeomAnalyzerI *> (_root_geom);
    
}

void MindGenieEventInterface::initGenieFluxDriver(){
  /* 
  int all_nu[4] = {12, -12, 14, -14};
  for (int i=0; i<sizeof(all_nu)/sizeof(int); i++)
    _fluxParticles.push_back(all_nu[i]);
  */
  _ntuple_flux_driver = new flux::GSimpleNtpFlux;
  _ntuple_flux_driver->LoadBeamSimData(_musr_flux_tree, 
				       MindConfigService::Instance().FileName());

  // _flux_driver = dynamic_cast<GFluxI*> (_ntuple_flux_driver)
}

void MindGenieEventInterface::initGenieMCjobDriver(){
  _mcj_driver = new GMCJDriver;
  _mcj_driver->UseFluxDriver(_ntuple_flux_driver);
  _mcj_driver->UseGeomAnalyzer(_geom_anal);
  _mcj_driver->Configure();
  _mcj_driver->UseSplines();
  // _mcj_driver->ForceSingleProbScale();
}

void MindGenieEventInterface::initGenieMCjobMonitor(){
  _mcjmonitor = new GMCJMonitor(_runNumber);
}

void MindGenieEventInterface::initGenieMessenger(){
  Messenger* msg = Messenger::Instance();
  msg->SetPriorityLevel("GHEP",pERROR);
  msg->SetPriorityLevel("xSecSplit",pERROR);
}

G4PrimaryParticle* 
MindGenieEventInterface::CreatePrimaryParticle(GHepParticle& part, 
					       G4int PDG, bool primLep)
{
  G4ParticleDefinition* particle_def;
  
  if ( PDG < 1000000000 )
    particle_def =
      G4ParticleTable::GetParticleTable()->FindParticle( PDG );
  else
    particle_def =
      G4ParticleTable::GetParticleTable()->GetIon( part.Z(), part.A(), 0 );
  
  G4double mass   = particle_def->GetPDGMass();
  
  G4PrimaryParticle* particle =  
    new G4PrimaryParticle(particle_def, part.Px()*GeV, part.Py()*GeV, part.Pz()*GeV);
  
  if ( primLep ){
    MindLookupTable::Instance().set_primLep( PDG, part.Px()*GeV, part.Py()*GeV, part.Pz()*GeV);
    /* cout<<"px = "<<part.Px()*GeV
	<<" py = "<<part.Py()*GeV
222	<<" pz = "<<part.Pz()*GeV<<endl; */
  }
  else {
    //Add to 'hadronic vector' GeV.
    _had4P[0] += part.Px();
    _had4P[1] += part.Py();
    _had4P[2] += part.Pz();
    _had4P[3] += part.E();
  } 
  /*
  cout<<"Hadron px = "<<_had4P[0]
      <<" py = "<<_had4P[1]
      <<" pz = "<<_had4P[2]<<endl;
  */
  particle->SetMass( mass );
  particle->SetCharge( part.Charge() );
  //This should maybe be updated to allow polarisation recorded by genie?
  //Generally only set by genie for lepton but could be relevant.
  if ( part.PolzIsSet() ){
    TVector3 plz;
    part.GetPolarization( plz );
    particle->SetPolarization( plz[0], plz[1], plz[2] );
  } else
    particle->SetPolarization(0.,0.,0.);
  
  return particle;
}

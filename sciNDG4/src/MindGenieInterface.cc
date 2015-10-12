#include "MindGenieInterface.h"

// #include <CLHEP/Units/GlobalSystemOfUnits.h>

// #include <bhep/bhep_svc.h>

MindGenieInterface::MindGenieInterface()
{
  Initialize();
}

MindGenieInterface::~MindGenieInterface()
{
  //Delete and close the input files.
  Finalize();
}

void MindGenieInterface::Initialize()
{
  
  const MindParamStore& config = MindConfigService::Instance().Generation();
  
  _mcrec[0] = new NtpMCEventRecord();
  _mcrec[1] = new NtpMCEventRecord();

  //_mcrec[0] = 0;
  //_mcrec[1] = 0;

  //if (config.PeekSParam("active_material_data")) {
  G4String active_mat_file  = config.GetSParam("active_material_data");
  _genieFiles[0] = new TFile( active_mat_file.data(), "read" );
  _genieData[0] = dynamic_cast <TTree *>( _genieFiles[0]->Get("gtree")  );
  _genieData[0]->SetBranchAddress("gmcrec", &_mcrec[0]);
  
  // Only get Iron target events if there is iron.
  if (config.PeekSParam("passive_material_data")) {
    G4String passive_mat_file = config.GetSParam("passive_material_data");
    _genieFiles[1] = new TFile( passive_mat_file.data(), "read" );
    _genieData[1] = dynamic_cast <TTree *>(_genieFiles[1]->Get("gtree")  );
    _genieData[1]->SetBranchAddress("gmcrec", &_mcrec[1]);
  }

  _vtx_location = config.GetSParam("vertex_location");

  if ( _vtx_location == "FIXED" ){
    std::vector<double> inVert = config.GetVParam("fvert");
    _fvec = G4ThreeVector( inVert[0], inVert[1], inVert[2] );
  }

  _fsl_select = config.PeekIParam("FSL_Select") ? config.GetIParam("FSL_Select") : 0;

  _evCount[0] = 0; _evCount[1] = 0;
  if ( MindConfigService::Instance().Geometry().PeekIParam("TASD") )
    _tasd = MindConfigService::Instance().Geometry().GetIParam("TASD");
  else _tasd = false;

  _had4P = bhep::vdouble(4);
  _fspdg = bhep::vdouble(6);
  for(int i=0; i<6; i++) _fspdg[i] = 0;
  
}

void MindGenieInterface::Finalize()
{

  delete _genieData[0];
  _genieFiles[0]->Close();
  delete _genieFiles[0];

  if ( _genieFiles[1] ) {
    delete _genieData[1];
    _genieFiles[1]->Close();
    delete _genieFiles[1];
  }

}

void MindGenieInterface::GeneratePrimaryVertex(G4Event* event)
{
  G4int region_code = SelectVertexRegion();
  //EOF protection.
  if ( _evCount[region_code] >= _genieData[region_code]->GetEntries() )
     G4Exception("MindGenieInterface::GeneratePrimaryVertex",
		"GenieInterface_S001", RunMustBeAborted,
		"ERROR.- End-Of-File reached!");

  // get bhep transient event from singleton
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();
  
  // std::cout<<"Reading Tree at "<<_evCount[region_code]<<std::endl;
  G4int catcher = _genieData[region_code]->GetEntry( _evCount[region_code] );

  _genieData[region_code]->Show(_evCount[region_code]);
  // start of the event

  G4double time = 0.;
  for (bhep::vdouble::iterator it1=_had4P.begin();it1 != _had4P.end();it1++)
    (*it1) = 0.;
    
  MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
    G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  
  G4String region_name;
  
  if (region_code == 0) region_name = "ACTIVE";
  else region_name = "PASSIVE";
  
  G4ThreeVector position = detConstr->GetDetectorGeometry()->
    GetVertex( region_name );

  if ( _vtx_location == "FIXED" )
    position = _fvec;
  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);
  
  //Set Vertex to bhep event.
  bevt.set_vertex( position.x(), position.y(), position.z() );

  //Get the event.
  EventRecord & gEvent = *(_mcrec[region_code]->event);
  
  //Event information.
  Interaction* gInt = gEvent.Summary();
  GHepParticle* fsl = gEvent.FinalStatePrimaryLepton();
  
  bevt.add_property( "IntType", gInt->ProcInfo().AsString() );

  bool charmProd = gInt->ExclTag().IsCharmEvent();
  bevt.add_property( "Charm", charmProd );
  if ( charmProd )
    bevt.add_property( "CharmHad", gInt->ExclTag().CharmHadronPdg() );

  G4double Q2 = gInt->Kine().Q2(true);
  if ( Q2 != -99999 )//Protects against cases where Q2 is not set.
    bevt.add_property( "Q2", gInt->Kine().Q2(true)*GeV );
  bevt.add_property( "EngTrans", (gEvent.Probe()->E() - fsl->E())*GeV );
  //

  //Loop over particles.
  TObjArrayIter iter(&gEvent);
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
	// G4cout<<part->E()*GeV<<G4endl;
      }
      else {
	//Electron scattering events have neutrino, electron and nucleus here.
	if ( pdg >= 1000000000 ) {//Nucleus
	  bevt.add_property("nucType",   pdg);
	  bevt.add_property("nucEnergy", part->E()*GeV);
	} else {
	  bevt.add_property("intpart",   pdg);
	  bevt.add_property("partEnergy", part->E()*GeV);
	  // G4cout<<pdg<<"\t"<<part->E()*GeV<<G4endl;
	}

      }
      
    } else if ( pStatus == 1 && pdg < 2000000000 ) { //stable final state.
      if(fslcount<6)_fspdg[fslcount] = pdg;
      fslcount++;
      if ( part == fsl ){
	// G4cout<<part->Pdg()<<"\t"<<part->E()<<G4endl;
	vertex->SetPrimary(CreatePrimaryParticle( *part, pdg, true ));
      }
      else if(_fsl_select==0){
	// G4cout<<part->Pdg()<<"\t"<<part->E()<<G4endl;
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
  
  _evCount[region_code]++;
  _mcrec[region_code]->Clear();

}

G4int MindGenieInterface::SelectVertexRegion()
{
  if      ( _tasd ) return 0;

  if      (_vtx_location == "ACTIVE") return 0;
  
  else if (_vtx_location == "PASSIVE") return 1;

  else if (_vtx_location == "RANDOM") {
    
    // Randomly select whether vertex should be located in
    // passive or active material
    
    static G4double passive_target_prob = GetTargetProb();
    
    if (G4UniformRand() <= passive_target_prob) return 1;
    else return 0;
  }

  else if (_vtx_location == "FIXED") {
    if ( _evCount[0] == 0 && _evCount[1] == 0 ){
    MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
      G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    _fvecReg = detConstr->
      GetDetectorGeometry()->GetRegion( _fvec[2] );
    }

    return _fvecReg;
  }
}

G4double MindGenieInterface::GetTargetProb()
{
  MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
    G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  
  return (detConstr->GetDetectorGeometry()->GetPassiveTargetProb());
}

G4PrimaryParticle* 
MindGenieInterface::CreatePrimaryParticle(GHepParticle& part, G4int PDG, bool primLep)
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

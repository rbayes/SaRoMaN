// ----------------------------------------------------------------------------
//  $Id: MindNuanceInterface.cc 543 2014-11-01 23:03:10Z  $
//
//  Authors: A Laing <a.laing@physics.gla.ac.uk>
//           J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 21 May 2009 
//
//  Copyright (c) 2009 - IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindNuanceInterface.h"

//#include <bhep/bhep_svc.h>



MindNuanceInterface::MindNuanceInterface()
{
  Initialize();
}



void MindNuanceInterface::Initialize()
{
  const MindParamStore& config = MindConfigService::Instance().Generation();

  G4String active_mat_file  = config.GetSParam("active_material_data");
  _nuanceFiles[0].open(active_mat_file.data());
  // Only get Iron target events if there is iron.
  if (config.PeekSParam("passive_material_data")) {
    G4String passive_mat_file = config.GetSParam("passive_material_data");
    _nuanceFiles[1].open(passive_mat_file.data());
  }
  
  _vtx_location = config.GetSParam("vertex_location");

  if ( _vtx_location == "FIXED" ){
    std::vector<double> inVert = config.GetVParam("fvert");
    _fvec = G4ThreeVector( inVert[0], inVert[1], inVert[2] );
  }

  if ( MindConfigService::Instance().Geometry().PeekIParam("TASD") )
    _tasd = (bool)MindConfigService::Instance().Geometry().GetIParam("TASD");
  else _tasd = false;

  _runStart = false;

  _had4P = bhep::vdouble(4);
}



void MindNuanceInterface::GeneratePrimaryVertex(G4Event* event)
{
  G4int region_code = SelectVertexRegion();

  //G4ThreeVector position(0.,0.,-18420);
  G4int NuPDG;
  
  // get bhep transient event from singleton
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();

  G4int track_type;
  G4String region, dollar, line_type;
  
  // read a line and check type
  _nuanceFiles[region_code] >> dollar >> line_type;
  
  if (line_type == "begin") {  // start of the event

    G4double time = 0.;
    
    MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
      G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    
    G4String region_name;

    if (region_code == 0) region_name = "ACTIVE"; 
    else region_name = "PASSIVE";

    //if (_vtx_location == "FIXED")
    //position = G4ThreeVector(0., 0., -18420);
    //else
    G4ThreeVector position = detConstr->GetDetectorGeometry()->
      GetVertex(region_name);

    if ( _vtx_location == "FIXED" )
      position = _fvec;
    
    G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);

    //Set Vertex to bhep event.
    bevt.set_vertex( position.x(), position.y(), position.z() );
    
    do {  // read and process lines until 'end' is reached  
      
      _nuanceFiles[region_code] >> dollar >> line_type;
      
      if (line_type == "nuance") {  // event interaction mode
	G4int interaction_code;
	_nuanceFiles[region_code] >> interaction_code;

	// store interaction type in bhep transient event as a property
	bevt.add_property("IntType", InteractionType(interaction_code));
	
	
      }
      
      else if (line_type == "track") {  // particle info ............
	
	G4int pdg;
	G4double energy, thetaX, thetaY, thetaZ;
	_nuanceFiles[region_code] >> pdg >> energy >> thetaX >> thetaY >> thetaZ >> track_type;
	
	if (track_type == 0) {  // final-state particle
	  G4ThreeVector direction(thetaX, thetaY, thetaZ);
	  if ( abs(pdg) == NuPDG || abs(pdg) == NuPDG-1 )
	    vertex->SetPrimary(CreatePrimaryParticle(pdg, energy, direction, true));
	  else
	    vertex->SetPrimary(CreatePrimaryParticle(pdg, energy, direction));
	}
	
	else if (track_type == -1) {  // initial-state particle

	  if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) {
	    NuPDG = pdg;
	    bevt.add_property("nuType",   pdg);
	    bevt.add_property("nuEnergy", energy*MeV);
	  }
	  else {
	    bevt.add_property("nucType",   pdg);
	    bevt.add_property("nucEnergy", energy*MeV);
	  }
	}
      } // end particle info ........................................

      else if (line_type == "info") {
	// obsolete/irrelevant infomation that we skip
	_nuanceFiles[region_code].ignore( 100, '\n' );
      }
      
      else if (line_type == "vertex") {
	// nuance event vertex that we don't use
	_nuanceFiles[region_code].ignore( 100, '\n' );
      }
    } while (line_type != "end");

    bevt.add_property("had4vec", _had4P);
    
    event->AddPrimaryVertex(vertex);
    
  }
  
  else if (line_type == "stop") {
    G4Exception("NuanceInterface::GeneratePrimaryVertex",
		"NuanceInterface_S001", RunMustBeAborted,
		"ERROR.- End-Of-File reached!");
  } 
  else {
    G4cout << "Bad line: " << line_type << G4endl;
    G4Exception("NuanceInterface::GeneratePrimaryVertex",
		"NuanceInterface_S002", RunMustBeAborted,
		"ERROR.- Unrecognized line type!");
  }
  
}



G4PrimaryParticle* 
MindNuanceInterface::CreatePrimaryParticle(G4int PDG, G4double eng, G4ThreeVector dir, bool pLep)
{
  G4double pmom, px, py, pz;
  G4int ionA, ionZ;
  //Ion protection, maybe needs fixed.
  G4ParticleDefinition* particle_def;
  
  if ( PDG < 6000 )
    particle_def =
      G4ParticleTable::GetParticleTable()->FindParticle(PDG);
  else {
    ionZ = PDG / 1000;
    ionA = PDG % 1000;
    particle_def =
      G4ParticleTable::GetParticleTable()->GetIon( ionZ, ionA, 0);
  }

  G4double mass   = particle_def->GetPDGMass();
  G4double charge = particle_def->GetPDGCharge();
  
  pmom = std::sqrt( eng*eng - mass*mass );
  px = pmom * dir.x();
  py = pmom * dir.y();
  pz = pmom * dir.z();
  
  G4PrimaryParticle* particle =  
    new G4PrimaryParticle(particle_def, px, py, pz);

  if ( pLep )
    MindLookupTable::Instance().set_primLep( PDG, px, py, pz);
  else {
    //Add to 'hadronic vector' GeV.
    _had4P[0] += px;
    _had4P[1] += py;
    _had4P[2] += pz;
    _had4P[3] += eng;
  } 
  
  particle->SetMass(mass);
  particle->SetCharge(charge);
  particle->SetPolarization(0.,0.,0.);

  return particle;
}



G4int MindNuanceInterface::SelectVertexRegion()
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
    if ( !_runStart ){
      MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
	G4RunManager::GetRunManager()->GetUserDetectorConstruction();
      _fvecReg = detConstr->
	GetDetectorGeometry()->GetRegion( _fvec.z() );
      _runStart = true; 
    }

    return _fvecReg;
  }
}



/// Convert nuance interaction code to a more user-friendly name
///
G4String MindNuanceInterface::InteractionType(G4int code)
{
  G4String name;
  
  if      ( abs(code) == 1 )  name = "CCQE"; 
  else if ( abs(code) == 2 )  name = "NCQE";
  else if ( abs(code) >= 3  && abs(code) <= 16 ) name = "1piRes";
  else if ( abs(code) >= 17 && abs(code) <= 90 ) name = "miscRes";
  else if ( abs(code) == 91 ) name = "CCDIS";
  else if ( abs(code) == 92 ) name = "NCDIS";
  else if ( abs(code) == 95 ) name = "CabQE";
  else if ( abs(code) == 96 ) name = "NCpi";
  else if ( abs(code) == 97 ) name = "CCpi";
  else if ( abs(code) == 98 ) name = "eEL";
  else if ( abs(code) == 99 ) name = "muINVe";
  else name = "unknown";

  return name;
}


G4double MindNuanceInterface::GetTargetProb()
{
  MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
    G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  
  return (detConstr->GetDetectorGeometry()->GetPassiveTargetProb());
}

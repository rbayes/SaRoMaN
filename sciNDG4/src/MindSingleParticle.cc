// ----------------------------------------------------------------------------
//  $Id: MindSingleParticle.cc 543 2014-11-01 23:03:10Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 17 March 2009
//  
//  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindSingleParticle.h"


MindSingleParticle::MindSingleParticle()
{
  Initialize();
}



void MindSingleParticle::Initialize()
{
  const MindParamStore& config = MindConfigService::Instance().Generation();

  // Set particle definition (i.e., type) by name
  G4String particle_name = config.GetSParam("particle_name");
  G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
  _particle_definition = particle_table->FindParticle(particle_name);

  // Set particle basic properties from definition
  _mass   = _particle_definition->GetPDGMass();
  _charge = _particle_definition->GetPDGCharge();

  _energy_min = config.GetDParam("energy_min") * GeV;
  _energy_max = config.GetDParam("energy_max") * GeV;

  if(config.PeekDParam("costh_min"))
    _costh_min = config.GetDParam("costh_min");
  else
    _costh_min = 1.0;
  if(config.PeekDParam("costh_max"))
    _costh_max = config.GetDParam("costh_max");
  else
    _costh_max = 1.0;

  _vtx_location = config.GetSParam("vertex_location");

  if ( _vtx_location == "FIXED" || _vtx_location == "GAUSS"){
    std::vector<double> inVert = config.GetVParam("fvert");
    _fvec = G4ThreeVector( inVert[0], inVert[1], inVert[2] );
  
    if (_vtx_location == "GAUSS"){
      _RMS = config.GetVParam("bspot");
    }
  }
  else {
    _fvec = G4ThreeVector( 0., 0., -2000.);
  }

  _evCount[0] = 0; _evCount[1] = 0;

  if ( MindConfigService::Instance().Geometry().PeekIParam("TASD") )
    _tasd = MindConfigService::Instance().Geometry().GetIParam("TASD");
  else _tasd = false;
  
  _had4P = bhep::vdouble(4);
  _fspdg = bhep::vdouble(6);
  for(int i=0; i<6; i++) _fspdg[i] = 0;
}



void MindSingleParticle::GeneratePrimaryVertex(G4Event* event)
{
  // G4int region_code = SelectVertexRegion();
  // particle is generated at time zero
  G4double time = 0.;
  
  for (bhep::vdouble::iterator it1=_had4P.begin();it1 != _had4P.end();it1++)
    (*it1) = 0.;

  // G4double length = 
  //   MindConfigService::Instance().Geometry().GetDParam("MIND_z") * mm;
  // G4ThreeVector particle_position = G4ThreeVector(0., 0., -length/2.);
  
    // get bhep transient event from singleton
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();

    
  MindDetectorConstruction* detConstr = (MindDetectorConstruction*) 
    G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  
  // G4String region_name;
  
  // if (region_code == 0) region_name = "ACTIVE";
  // else region_name = "PASSIVE";
  
  G4ThreeVector position;
   
  if( _vtx_location == "GAUSS"){
    double x = G4RandGauss::shoot(_fvec[0], _RMS[0]);
    double y = G4RandGauss::shoot(_fvec[1], _RMS[1]);
    position = _fvec;
    position[0] = x;
    position[1] = y;
  }
  // G4cout<<position.x()<<"\t"<<position.y()<<"\t"<<position.z()<<"\n";
  if ( _vtx_location == "FIXED" )
    position = _fvec;
  
  // position.setZ(-length/2.);

  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);
  
  //Set Vertex to bhep event.
  bevt.set_vertex( position.x(), position.y(), position.z() );
  // Added to allow the analysis
  // Not really true, but sufficient to make the trains run
  bevt.add_property( "IntType", "CCQE" );
  bevt.add_property( "Charm", false );
  // generate random energy in [E_min, E_max]
  G4double kinetic_energy = GenerateRandomEnergy(_energy_min, _energy_max);

  // forward direction
  G4ThreeVector momentum_direction(0., 0., 1.);
  if(_costh_min != 1.)
    momentum_direction = GenerateRandomDirection();
  
  // set cartesian components of momentum from kinetic energy
  G4double energy = kinetic_energy + _mass; // total energy
  G4double pmom = std::sqrt(energy*energy - _mass*_mass); // momentum module
  G4double px = pmom * momentum_direction.x();
  G4double py = pmom * momentum_direction.y();
  G4double pz = pmom * momentum_direction.z();

  // create new particle and add it to the primary vertex
  G4PrimaryParticle* particle = 
    new G4PrimaryParticle(_particle_definition, px, py, pz);
  particle->SetMass(_mass);
  particle->SetCharge(_charge);
  particle->SetPolarization(0.,0.,0.);
  vertex->SetPrimary(particle);
  
  bevt.add_property("partEnergy", energy);
  bevt.add_property("intpart",_particle_definition->GetPDGEncoding());
  // The following lines are false but it is a convienient fiction
  bevt.add_property("nuEnergy", energy);
  bevt.add_property("nupart",_particle_definition->GetPDGEncoding());
  // add primary vertex to event
  
  _fspdg[0] = _particle_definition->GetPDGEncoding();
  
  bevt.add_property("had4vec", _had4P);
  bevt.add_property("npip",0);
  bevt.add_property("npi0",0);
  bevt.add_property("npim",0);
  bevt.add_property("nkp",0);
  bevt.add_property("nk0",0);
  bevt.add_property("nkm",0);
  bevt.add_property("hadEInNucleus", 0.0);
  bevt.add_property("fspdg", _fspdg);

  event->AddPrimaryVertex(vertex);
}


G4ThreeVector MindSingleParticle::GenerateRandomDirection()
{
  // Random distribution uniform in solid angle
  G4double cos_theta = _costh_min == _costh_max ? _costh_min: 
    (_costh_max - _costh_min) * G4UniformRand() - _costh_min;
  G4double sin_theta = std::sqrt(1. - cos_theta*cos_theta);
  G4double phi = twopi * G4UniformRand();

  G4double ux = sin_theta * std::cos(phi);
  G4double uy = sin_theta * std::sin(phi);
  G4double uz = cos_theta;
  
  return G4ThreeVector(ux, uy, uz);
}



G4double MindSingleParticle::GenerateRandomEnergy(G4double min, G4double max)
{
  if(_vtx_location == "GAUSS"){
  
    if(abs(_particle_definition->GetPDGEncoding())==211){
      double pmin = 1 - exp(-min/2. / GeV);
      double pmax = 1 - exp(-max/2. / GeV);
      return -2. * GeV * log (1 - pmin + G4UniformRand() * (pmax - pmin));
    }
    if(abs(_particle_definition->GetPDGEncoding())==13){
      double pmin = exp(-min/4./ GeV);
      double pmax = exp(-max/4./ GeV);
      return -4. * GeV * log (pmin + G4UniformRand() * (pmax - pmin));
    }
    else return (G4UniformRand() * (max - min) + min);
  }
  else  return (G4UniformRand() * (max - min) + min);
}


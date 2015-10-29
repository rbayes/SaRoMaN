// ----------------------------------------------------------------------------
//  $Id: SciNearDetectorGeometry.cc 467 2011-11-25 16:00:27Z ryan $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 14 Jun 2008
//  
//  Copyright (c) 2008, 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "SciNearDetectorGeometry.h"
#include "MindConfigService.h"
#include "MindParamStore.h"
#include "MindSD.h"
#include "MindField.hh"
#include "MindFieldMapR.hh"
#include "babyMindField.hh"

#include <globals.hh>
#include <G4Box.hh>
#include <G4ExtrudedSolid.hh>
#include <G4Tubs.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4PVReplica.hh>
#include <Randomize.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4TransportationManager.hh>
#include <G4SDManager.hh>
#include <G4ChordFinder.hh>
#include <G4UserLimits.hh>

//#include <G4PVParameterised.hh>

#include <bhep/utilities.h>


SciNearDetectorGeometry::SciNearDetectorGeometry()
{
  SetInputParameters();
  CalculatePassiveTargetProb();
  _logVolume = DefineDetector();
}

void SciNearDetectorGeometry::SetInputParameters()
{
  const MindParamStore& config =
    MindConfigService::Instance().Geometry();
  
  // Parameters for MIND muon catcher.
  _piece_width       = config.GetDParam("MIND_x") * m;
  _piece_height      = config.GetDParam("MIND_y") * m;
  _mind_length       = config.GetDParam("MIND_z") * m;
  _ear_width         = config.PeekDParam("ear_width")?
    config.GetDParam("ear_width") *m:0.1*mm;
  _ear_height        = config.PeekDParam("ear_height")?
    config.GetDParam("ear_height") *m:_piece_height/2.;
  _bore_diameter     = config.PeekDParam("bore_diameter")?
    config.GetDParam("bore_diameter")*m: 100*mm;

  // Parameters for a vertex detector.
  _vertex_width      = config.GetDParam("vertex_x") * m;
  _vertex_height     = config.GetDParam("vertex_y") * m;
  _vertex_depth      = config.GetDParam("vertex_z") * m;
  // _cal_depth         = config.GetDParam("cdepth") * mm;

  if(_vertex_width == 0 || _vertex_height == 0 || _vertex_depth == 0)
    {
      _isVertexDet = 0;
    }
  else
    {
      _isVertexDet = 1;
    }

  // Use the same active thickness for the vertex detector and calorimeter
  _active_thickness  = config.GetDParam("active_thickness") * cm;
  _passive_thickness = config.GetDParam("passive_thickness") * cm;
  _bracing_thickness = config.PeekDParam("bracing_thickness") * cm ?
    config.GetDParam("bracing_thickness") * cm : 0.0 * mm;

  _number_active = config.GetIParam("active_layers");

  IsOctagonal = config.PeekIParam("IsOctagonal")?
    config.GetIParam("IsOctagonal"): 1;
  //if(IsOctagonal) _bracing_thickness = 0.0 * cm;
  
  _piece_length = 0;
  _TASDm_length = 0;
  G4String gap_num;
  for (G4int ilayer = 1;ilayer <= _number_active+1;ilayer++){
    gap_num = bhep::to_string(ilayer);

    if ( config.PeekDParam("gap"+gap_num) ){
      _gaps.push_back( config.GetDParam("gap"+gap_num) * mm );
      _piece_length += _gaps[ilayer-1];
      if( ilayer < _number_active+1 ) 
	// TASD has no passive layers, so one less gap
	_TASDm_length += _gaps[ilayer-1];
    } else
      G4Exception("SciNearDetectorGeometry::SetInputParameters",
		  "Geometry_S001",RunMustBeAborted,
		  "ERROR.- Missing gap size!");
  }
  
  _piece_length += _passive_thickness + _number_active*_active_thickness 
    + _number_active * _bracing_thickness;

  _TASDm_length += _number_active*_active_thickness;

  _active_displacement = std::vector<G4double> (_number_active, 0);
  _brace_displacement = std::vector<G4double> (_number_active, 0);
  _detector_length = _mind_length + _vertex_depth;
  _number_pieces = (G4int)ceil(_mind_length / _piece_length);
  _number_TASDpieces  = (G4int)ceil((_vertex_depth) / _TASDm_length);
  _detector_length = _number_pieces * _piece_length + _number_TASDpieces * _TASDm_length;
  _mind_length = _number_pieces * _piece_length;
  _vertex_depth = _number_TASDpieces * _TASDm_length;

  _active_material  = config.GetSParam("active_material");
  _passive_material = config.GetSParam("passive_material");
  _bracing_material = config.PeekSParam("bracing_material") ?
    config.GetSParam("bracing_material") : "G4_Al";

  //Are we making a tasd? here we always assume yes 
  // but not in the intended way.
  if ( config.PeekIParam("TASD") )
    _tasd = config.GetIParam("TASD");
  else
    _tasd = 0;
  
  // And again for this detector we know the fields.
  /*
  _isUniform = false;
  if ( config.PeekIParam("isUniform") ){
    int ubit = config.GetIParam("isUniform");
    if ( ubit == 1 ) _isUniform = true;
  }
  */
}

G4LogicalVolume* SciNearDetectorGeometry::DefineDetector()
{
  // PIECE ..........................................................
  // It is the basic element of the calorimeter. It consists of a
  // PASSIVE absorber layer and two ACTIVE detection layers.
  // It is filled with air.
  
  const int nz = 2;
  double * z = new double[nz];
  double * rin = new double[nz];
  double * rout = new double[nz];
  const double Pi = 4*atan(1);
  double lh = _piece_height/2. * tan(Pi/8);
  double lv = _piece_width/2. * tan(Pi/8);
  double llv = _piece_width/2. + _ear_width;
  double llh = lv - _ear_width;
  double hst = _ear_height - _piece_height/2.;
  
  // 


  std::vector<G4TwoVector> MindSection;
  
  if(IsOctagonal == 1){ 
    // Then use the parameters to define the eared octagon
    MindSection.push_back( G4TwoVector(lh, _piece_height/2.) );
    MindSection.push_back( G4TwoVector(llv, llh) );
    MindSection.push_back( G4TwoVector(llv, hst) );
    MindSection.push_back( G4TwoVector(_piece_width/2., hst) );
    MindSection.push_back( G4TwoVector(_piece_width/2., -lv) );
    MindSection.push_back( G4TwoVector(lh, -_piece_height/2.) );
    MindSection.push_back( G4TwoVector(-lh, -_piece_height/2.) );
    MindSection.push_back( G4TwoVector(-_piece_width/2., -lv) );
    MindSection.push_back( G4TwoVector(-_piece_width/2., hst) );
    MindSection.push_back( G4TwoVector(-llv, hst) );
    MindSection.push_back( G4TwoVector(-llv, llh) );
    MindSection.push_back( G4TwoVector(-lh, _piece_height/2.) );
  }
  else if(IsOctagonal == 2){
    // Then use the parameters to define the doubled octogon for LBNO
    // No Ears for now. 
    MindSection.push_back( G4TwoVector(lh + _piece_width/4., _piece_height/2.) );
    MindSection.push_back( G4TwoVector(_piece_width/2., lh) );
    MindSection.push_back( G4TwoVector(_piece_width/2., -lh) );
    MindSection.push_back( G4TwoVector(lh + _piece_width/4., -_piece_height/2.) );
    MindSection.push_back( G4TwoVector(-lh - _piece_width/4., -_piece_height/2.) );
    MindSection.push_back( G4TwoVector(-_piece_width/2., -lh) );
    MindSection.push_back( G4TwoVector(-_piece_width/2.,  lh) );
    MindSection.push_back( G4TwoVector(-lh-_piece_width/4., _piece_height/2.) );
  }
  /*
  else { // use a cylindrical cross-section
    // verticle and horizontal directions should be the same, but
    // assuming that they are not, generalize to an ellipse.
    for(double th=0; th < 2*pi; th+=pi/90.){
      double x_t = _piece_width * cos(th);
      double y_t = _piece_height * sqrt( 1 - cos(th)*cos(th));
      if(th < pi)
	MindSection.push_back( G4TwoVector(x_t, y_t) );
      else
	MindSection.push_back( G4TwoVector(x_t,-y_t) );
    }
    }*/
  G4LogicalVolume* piece_logic;
  
  if(IsOctagonal==1 || IsOctagonal==2){
    G4ExtrudedSolid* piece_solid = 
      new G4ExtrudedSolid("PIECE",MindSection,_piece_length/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
    piece_logic = 
      new G4LogicalVolume(piece_solid, G4Material::GetMaterial("G4_AIR"), 
			  "PIECE", 0, 0, 0, true);
  }
  else if(IsOctagonal == 3) {
    G4Box* piece_solid = 
      new G4Box("PIECE", (_piece_width+_ear_width)/2., (_piece_height+_ear_height)/2., _piece_length/2.);
    piece_logic = 
      new G4LogicalVolume(piece_solid, G4Material::GetMaterial("G4_AIR"), 
			  "PIECE", 0, 0, 0, true);
  } 
  else {
    G4Tubs* piece_solid = 
      new G4Tubs("PIECE", 0.*cm, _piece_height/2., _piece_length/2., 0.0*deg, 360.*deg);
    piece_logic = 
      new G4LogicalVolume(piece_solid, G4Material::GetMaterial("G4_AIR"), 
			  "PIECE", 0, 0, 0, true);
  }
  /*
  z[0] = 0.;
  z[1] = _piece_length/2.;
  rin[0] = 0.0;
  rin[1] = rin[0];
  rout[0] = 
  rout[1] = rout[0];

  G4Polyhedra* piece_solid = 
    new G4Polyhedra("PIECE", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
  */
  
  
  // PASSIVE ........................................................
 
  G4LogicalVolume* passive_logic;

  if(IsOctagonal==1){
    G4ExtrudedSolid* FePlate_solid = 
      new G4ExtrudedSolid("PASSIVE",MindSection,_passive_thickness/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
  
    G4Tubs* STL_holeFe = new G4Tubs("STL_holeFe", 0.*cm, _bore_diameter/2., _passive_thickness/2., 0.0*deg, 360.*deg);
    
    G4SubtractionSolid* passive_solid = new G4SubtractionSolid("PASSIVE",FePlate_solid,STL_holeFe,0, G4ThreeVector());
  // z[1] = _passive_thickness;
  // G4Polyhedra* passive_solid = 
  //   new G4Polyhedra("PASSIVE", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
    
    passive_logic = new G4LogicalVolume(passive_solid, G4Material::GetMaterial(_passive_material), 
					"PASSIVE", 0, 0, 0, true);
  }
  else if(IsOctagonal==2){
    G4ExtrudedSolid* FePlate_solid = 
      new G4ExtrudedSolid("PASSIVE",MindSection,_passive_thickness/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
  
    G4Tubs* STL_holeFe = new G4Tubs("STL_holeFe", 0.*cm, _bore_diameter/2., _passive_thickness/2., 0.0*deg, 360.*deg);
    
    G4SubtractionSolid* temppas_solid = new G4SubtractionSolid("PASSIVE",FePlate_solid,STL_holeFe,0, G4ThreeVector(_piece_width/4.,0.,0.)); 
    G4SubtractionSolid* passive_solid = new G4SubtractionSolid("PASSIVE",temppas_solid,STL_holeFe,0, G4ThreeVector(-_piece_width/4.,0.,0.));
  // z[1] = _passive_thickness;
  // G4Polyhedra* passive_solid = 
  //   new G4Polyhedra("PASSIVE", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
    
    passive_logic = new G4LogicalVolume(passive_solid, G4Material::GetMaterial(_passive_material), 
					"PASSIVE", 0, 0, 0, true);
  }
  else if(IsOctagonal==3){
    G4Box* passive_solid = 
      new G4Box("PASSIVE", (_piece_width+_ear_width)/2., (_piece_height+_ear_height)/2., _passive_thickness/2.);
    
    passive_logic = new G4LogicalVolume(passive_solid, G4Material::GetMaterial(_passive_material), 
					"PASSIVE", 0, 0, 0, true);
  }
  else{
    G4Tubs* passive_solid = 
      new G4Tubs("PASSIVE", _bore_diameter/2., _piece_height/2., _passive_thickness/2., 0.0*deg, 360.*deg);
    
    passive_logic = new G4LogicalVolume(passive_solid, G4Material::GetMaterial(_passive_material), 
					"PASSIVE", 0, 0, 0, true);
  }

  //G4LogicalVolume* Fe_logic = 
  //  new G4LogicalVolume(FePlate_solid, G4Material::GetMaterial(_passive_material),
  //			"PASSIVE", 0, 0, 0, true);

  //SetMinimum Kinetic energy for tracking. from parameter file.
  const MindParamStore& config = MindConfigService::Instance().Physics();
  G4double minKin = config.GetDParam("minimum_kinEng") * MeV;

  passive_logic->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX, minKin));
  
  _passive_displacement = -_piece_length/2. + _passive_thickness/2.;
  
  G4PVPlacement* passive_physi = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, _passive_displacement), 
		      passive_logic, "PASSIVE", piece_logic, false, 0);
  
  // SetMagneticField( *passive_logic );
  // ACTIVE ........................................................
  
  G4LogicalVolume* active_logic;
  G4LogicalVolume* bracing_logic;
  if(IsOctagonal==1){
    G4ExtrudedSolid* SciPlate_solid = 
      new G4ExtrudedSolid("SciPlate_solid",MindSection,_active_thickness/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
    
    G4Tubs* STL_holeSci = new G4Tubs("STL_holeSci", 0.*cm, _bore_diameter/2., _active_thickness/2., 0.0*deg, 360.*deg);
    
    G4SubtractionSolid* active_solid = new G4SubtractionSolid("ACTIVE",SciPlate_solid,STL_holeSci,0, G4ThreeVector());
    
    //z[1] = _active_thickness;
    //G4Polyhedra* active_solid = 
    //  new G4Polyhedra("ACTIVE", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
    
    active_logic = 
      new G4LogicalVolume(active_solid, G4Material::GetMaterial(_active_material), 
			  "ACTIVE", 0, 0, 0, true);
  }
  else if(IsOctagonal==2){
    G4ExtrudedSolid* SciPlate_solid = 
      new G4ExtrudedSolid("SciPlate_solid",MindSection,_active_thickness/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
    
    G4Tubs* STL_holeSci = new G4Tubs("STL_holeSci", 0.*cm, _bore_diameter/2., _active_thickness/2., 0.0*deg, 360.*deg);
    
    G4SubtractionSolid* temp_solid = new G4SubtractionSolid("ACTIVE",SciPlate_solid,STL_holeSci,0, G4ThreeVector(_piece_width/4.,0.,0.));
    G4SubtractionSolid* active_solid = new G4SubtractionSolid("ACTIVE",temp_solid,STL_holeSci,0, G4ThreeVector(-_piece_width/4.,0.,0.));
    
    //z[1] = _active_thickness;
    //G4Polyhedra* active_solid = 
    //  new G4Polyhedra("ACTIVE", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
    
    active_logic = 
      new G4LogicalVolume(active_solid, G4Material::GetMaterial(_active_material), 
			  "ACTIVE", 0, 0, 0, true);
  }
  else if(IsOctagonal==3){
    G4Box* active_solid = 
      new G4Box("SciPlate_solid", _piece_width/2., _piece_height/2., _active_thickness/2.); 
    
    active_logic = 
      new G4LogicalVolume(active_solid, G4Material::GetMaterial(_active_material), 
			  "ACTIVE", 0, 0, 0, true);
    if(_bracing_thickness > 0.0){
      G4Box* bracing_solid = 
	new G4Box("Bracing_solid",_piece_width/2., _piece_height/2., _bracing_thickness/2.); 
      bracing_logic = 
	new G4LogicalVolume(bracing_solid, G4Material::GetMaterial(_bracing_material), 
			    "BRACE", 0, 0, 0, true);
    }

  }
  else {
    G4Tubs* active_solid = 
      new G4Tubs("SciPlate_solid", _bore_diameter/2., _piece_height/2., _active_thickness/2., 0.*deg, 360.*deg); 
    
    active_logic = 
      new G4LogicalVolume(active_solid, G4Material::GetMaterial(_active_material), 
			  "ACTIVE", 0, 0, 0, true);
    
    if(_bracing_thickness > 0.0){
      G4Tubs* bracing_solid = 
	new G4Tubs("Bracing_solid",_bore_diameter/2., _piece_height/2., _bracing_thickness/2., 0.*deg, 360.*deg); 
      bracing_logic = 
	new G4LogicalVolume(bracing_solid, G4Material::GetMaterial(_bracing_material), 
			    "BRACE", 0, 0, 0, true);
    }
    
  }
  //G4LogicalVolume* SciPlate_solid = 
  //  new G4LogicalVolume(SciPlate_solid, G4Material::GetMaterial(_active_material), 
  //			"ACTIVE", 0, 0, 0, true);

  //SetMinimum Kinetic energy for tracking.
  active_logic->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX, minKin));
  
  G4PVPlacement* active_physi;
  G4PVPlacement* brace_physi;

  G4String act_num;
  for (G4int iactive = 0;iactive < _number_active;iactive++){
    act_num = bhep::to_string( iactive+1 );

    _active_displacement[iactive] = -_piece_length/2. + _passive_thickness + _gaps[0] 
      + _active_thickness/2.;

    if ( iactive != 0 ){
      _active_displacement[iactive] += iactive * (_active_thickness + _bracing_thickness);
      for (G4int igap = 1;igap <= iactive;igap++)
	_active_displacement[iactive] += _gaps[igap];
    }

    active_physi = new G4PVPlacement(0, G4ThreeVector(0, 0, _active_displacement[iactive]), 
				   active_logic, "ACTIVE"+act_num, piece_logic, false, 0);
    if ( _bracing_thickness > 0.0 ){
      _brace_displacement[iactive] = -_piece_length/2. + _passive_thickness + _gaps[0] 
      + _active_thickness + _bracing_thickness/2. + _gaps[1]/2.;
      brace_physi = new G4PVPlacement(0, G4ThreeVector(0, 0, _brace_displacement[iactive]),
				      bracing_logic, "BRACE"+act_num, piece_logic, false, 0);
    }
  }
  
  // ACTIVE layers defined as sensitive detectors
  G4SDManager* SDMgr = G4SDManager::GetSDMpointer();
  G4String SD_name = "/mind/ACTIVE/";
  MindSD* MSD = new MindSD(SD_name);
  SDMgr->AddNewDetector(MSD);
  active_logic->SetSensitiveDetector(MSD);
  if ( _tasd == 1 ){
    G4cout << "Setting 'passive' as sensitive" << G4endl;
    passive_logic->SetSensitiveDetector(MSD);
  }

  // MIND .......................................................
  // It is a vacuum-filled volume that contains PIECES MIND. 
  G4LogicalVolume* mind_logic;
  if(IsOctagonal==1 || IsOctagonal==2){
    G4ExtrudedSolid* mind_solid = 
      new G4ExtrudedSolid("MIND",MindSection,_mind_length/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
    
    //z[1] = _mind_length;
    //G4Polyhedra* mind_solid = 
    //  new G4Polyhedra("MIND", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
    
    mind_logic =
      new G4LogicalVolume(mind_solid, G4Material::GetMaterial("G4_AIR"), 
			  "ACTIVE", 0, 0, 0, true);
  }
  else if(IsOctagonal==3){
    G4Box* mind_solid = 
      new G4Box("MIND", _piece_width/2., _piece_height/2., _mind_length/2.);
    
    mind_logic =
      new G4LogicalVolume(mind_solid, G4Material::GetMaterial("G4_AIR"), 
			  "ACTIVE", 0, 0, 0, true);
  }
  else {
    G4Tubs* mind_solid = 
      new G4Tubs("MIND",0.*cm, _piece_height/2., _mind_length/2., 0.*deg, 360*deg);
    
    mind_logic =
      new G4LogicalVolume(mind_solid, G4Material::GetMaterial("G4_AIR"), 
			  "ACTIVE", 0, 0, 0, true);
  }
  
  G4PVReplica* piece_physi = 
    new G4PVReplica("PIECE", piece_logic, mind_logic, 
		    kZAxis, _number_pieces, _piece_length);


  G4LogicalVolume* vertex_logic;

  if(_isVertexDet)
    {
      // Vertex detector layers
      G4Box* vlayer_solid = new G4Box("vlayer_solid", _vertex_width/2., _vertex_height/2., _active_thickness/2.);
      G4LogicalVolume* vlayer_logic = new G4LogicalVolume(vlayer_solid, 
							  G4Material::GetMaterial(_active_material),
							  "ACTIVE", 0, 0, 0, true);
      vlayer_logic->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX, minKin));
      vlayer_logic->SetSensitiveDetector(MSD);
      G4PVPlacement* vlayer_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0),vlayer_logic,"VLAYER", false, 0, true);
  
      // Now create the vertex detector
      G4Box* vertex_solid = new G4Box("vertex_solid", _vertex_width/2., _vertex_height/2., _vertex_depth/2.);
      vertex_logic = new G4LogicalVolume(vertex_solid, G4Material::GetMaterial("G4_AIR"),
							  "ACTIVE", 0, 0, 0,true);
      G4PVReplica* vlayer_physi = new G4PVReplica("LAYER", vlayer_logic, vertex_logic, kZAxis, 
						  _number_TASDpieces, _TASDm_length);
    }


  // DETECTOR .......................................................
  // It is a vacuum-filled volume that contains all PIECES. 
  G4LogicalVolume* detector_logic;
  if(IsOctagonal==1 || IsOctagonal==2){
    G4ExtrudedSolid* detector_solid = 
      new G4ExtrudedSolid("DETECTOR",MindSection,_detector_length/2.,
			  G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0);
    
    //z[1] = _detector_length;
    //G4Polyhedra* detector_solid = 
    //  new G4Polyhedra("DETECTOR", Pi/8., 2.0*Pi, 8, nz, z, rin, rout);
    
    detector_logic =
      new G4LogicalVolume(detector_solid, G4Material::GetMaterial("G4_AIR"), 
			  "ACTIVE", 0, 0, 0, true);
  }
  else if(IsOctagonal==3){
    G4Box* detector_solid = 
      new G4Box("DETECTOR", (_piece_width+_ear_width)/2., (_piece_height+_ear_height)/2., _detector_length/2.);
    
    detector_logic =
      new G4LogicalVolume(detector_solid, G4Material::GetMaterial("G4_AIR"), 
			  "ACTIVE", 0, 0, 0, true);
  }
  else {
    G4Tubs* detector_solid = 
      new G4Tubs("DETECTOR",0.*cm, _piece_height/2., _detector_length/2., 0.*deg, 360*deg);
    
    detector_logic =
      new G4LogicalVolume(detector_solid, G4Material::GetMaterial("G4_AIR"), 
			  "ACTIVE", 0, 0, 0, true);
  }
  G4VPhysicalVolume* mind_physi = new G4PVPlacement(0, G4ThreeVector(0, 0, _vertex_depth/2.),
						    mind_logic, "MIND", detector_logic, 0, 0);
  if(_isVertexDet)
    {
  G4VPhysicalVolume* vertex_physi = new G4PVPlacement(0, G4ThreeVector(0, 0, -_mind_length/2.),
						    vertex_logic, "VERTEX", detector_logic, 0, 0);
    }
  
  
    // The transmission line itself. Model it with solid copper for scattering
  /*
    G4Tubs* STL = new G4Tubs("STL", 0.*cm, 3.5*cm, _detector_length/2., 0.0*deg, 360.*deg);
    
    G4LogicalVolume* STL_logic = new G4LogicalVolume(STL, G4Material::GetMaterial("G4_Cu"),
						     "STL", 0, 0, 0, true);
    G4VPhysicalVolume* STL_physi = 
      new G4PVPlacement(0, G4ThreeVector(0,0,0), STL_logic, "STL", 0, false, 0);
  */
    

  // Visualization attributes ......................................
  
  G4VisAttributes* red  = new G4VisAttributes(G4Colour(1., 0., 0.));
  G4VisAttributes* blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  G4VisAttributes* green = new G4VisAttributes(G4Colour(0., 0.75, 0.25));
  
  //mind_logic->SetVisAttributes(red);
  mind_logic   ->SetVisAttributes(G4VisAttributes::Invisible);
  detector_logic   ->SetVisAttributes(G4VisAttributes::Invisible);
  piece_logic   ->SetVisAttributes(G4VisAttributes::Invisible);
  passive_logic   ->SetVisAttributes(G4VisAttributes::Invisible);
  passive_logic ->SetVisAttributes(red);
  active_logic  ->SetVisAttributes(green);
  if(_isVertexDet)
    {
  vertex_logic  ->SetVisAttributes(green);
    }
  // STL_logic  ->SetVisAttributes(green);
  
  if(_isUniform)
    SetMagneticField( *detector_logic );
  else {
    SetNullField( *detector_logic );
    SetMagneticField( *passive_logic );
    SetNullField( *active_logic );
    if(_isVertexDet)
      {
	SetDipoleField( *vertex_logic );
      }  
  }
  // gdm->Write("MIND_detector",detector_logic);

  return detector_logic;
}



void SciNearDetectorGeometry::SetNullField(G4LogicalVolume& detector_logic)
{
  // apply a global uniform magnetic field along Y axis
  G4FieldManager* fieldMgr = new G4FieldManager();
  _magField = 0;
  fieldMgr->SetDetectorField(_magField);
  fieldMgr->CreateChordFinder(_magField);
  fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);
  
  detector_logic.SetFieldManager( fieldMgr, false );
}


void SciNearDetectorGeometry::SetDipoleField(G4LogicalVolume& detector_logic)
{
  G4FieldManager* fieldMgr = new G4FieldManager();
  //if(_magField) delete magField; //delete the existing magn field
  const MindParamStore& config =
    MindConfigService::Instance().Geometry();

  std::cout<<"apply a regional uniform magnetic field along Y axis"<<std::endl;
  std::vector<G4double> fieldValue;  
  if ( config.PeekVParam("field") )
    fieldValue = config.GetVParam("field");
  
  double fieldScaling = +1.0;
  if(config.PeekDParam("FieldScaling"))
    fieldScaling = config.GetDParam("FieldScaling");
  
  G4UniformMagField* magField = 
    new G4UniformMagField(G4ThreeVector(fieldScaling*fieldValue[0]*tesla, 
					fieldScaling*fieldValue[1]*tesla, 
					fieldScaling*fieldValue[2]*tesla));
  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);
  fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);
  
  detector_logic.SetFieldManager( fieldMgr, true );
}

void SciNearDetectorGeometry::SetMagneticField(G4LogicalVolume& detector_logic)
{
  // apply a global uniform magnetic field along Y axis
  G4FieldManager* fieldMgr = new G4FieldManager();

  //if(_magField) delete magField; //delete the existing magn field
  const MindParamStore& config =
    MindConfigService::Instance().Geometry();

  std::vector<G4double> fieldValue;  
  if ( config.PeekVParam("field") )
    fieldValue = config.GetVParam("field");
  
  double fieldScaling = +1.0;
  if(config.PeekDParam("FieldScaling"))
    fieldScaling = config.GetDParam("FieldScaling");
  if (fieldValue.size() != 0) { 
    if(config.PeekSParam("FieldMap")){
      // A field map has been provided.  Should add capability to
      // scale field -- but would such a scaling be physical???  As it
      // is now the "field" parameter is just a switch to turn the
      // field on and off.
      G4String Bmap = config.GetSParam("FieldMap");
      // Declaration of the magnetic field map object
      MindFieldMapR* magField = new MindFieldMapR(Bmap,fieldScaling,
						  _passive_thickness, 
						  _number_active, 
						  _active_thickness);
      // Now to embed the field in the detector geometry
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
      fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);
    }
    else { // Field is not uniform. Use a toroidal field.
      std::cout<<"Using toroidal magnetic field"<<std::endl;
      double fieldScaling = +1.0;
      if(config.PeekDParam("FieldScaling"))
	fieldScaling = config.GetDParam("FieldScaling");
      bool isLBNO = IsOctagonal == 2 ? true:false;

      MindField* magField = new MindField(fieldScaling, _piece_height/2., _piece_width/2., isLBNO);
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
      fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);

    }
  }   
  else {
    G4cout << "No magnetic field" << G4endl;
    _magField = 0;
    fieldMgr->SetDetectorField(_magField);
  }

  detector_logic.SetFieldManager( fieldMgr, true );

}



G4ThreeVector SciNearDetectorGeometry::GetVertex(G4String region)
{
  // A calorimeter piece is selected randomly
  G4int piece_no = 0;

  if ( _tasd == 1 ){
    if ( G4UniformRand() < 1/(_number_active+1) )
      region = "PASSIVE";
    else region = "ACTIVE";
  }
  if ( _tasd == 2 ){
    region = "ACTIVE";
  }

  if (region == "PASSIVE") {
    piece_no = (G4int)ceil(_number_pieces * G4UniformRand());
    // Point in PASSIVE system of reference
    G4ThreeVector vtx = 
      RandomPointInBox(_piece_width, _piece_height, _passive_thickness, false);
    
    // Transformation to global system of reference
    // Displacement only in Z-axis direction; X,Y coordinates are ok.
    G4double Z = vtx.z() + _passive_displacement + 
      piece_no * _piece_length - _piece_length/2. - _detector_length/2. + _vertex_depth;
    
    vtx.setZ(Z);
    
    return vtx;
  }
  else if (region == "ACTIVE") {
    G4double vert_vol = _vertex_width * _vertex_height * _vertex_depth;
    G4double mind_scivol = _active_thickness*_vertex_height*_vertex_height * atan(1);
    G4double vert_prob = vert_vol / (vert_vol + mind_scivol);
    bool inVertexDet = G4UniformRand() < vert_prob || _tasd==2;

    piece_no = (G4int)ceil((inVertexDet?_number_TASDpieces:_number_pieces) 
			   * G4UniformRand());
    // Point in ACTIVE system of reference
    G4ThreeVector vtx = inVertexDet ?
      RandomPointInBox(_vertex_width, _vertex_height, _active_thickness, true):
      RandomPointInBox(_piece_width, _piece_height, _active_thickness, false);
    
    // Transformation to global system of reference
    // Displacement only in Z-axis direction; X,Y coordinates are ok.
    
    // Which ACTIVE. 1 or ACTIVE2? Selected randomly.
    G4int index = (G4int)floor( _number_active * G4UniformRand() );
    G4double displacement = _active_displacement[index];

    G4double Z = vtx.z() + displacement + 
      (inVertexDet ?
       piece_no * _TASDm_length - _TASDm_length/2. - _detector_length/2.:
       piece_no * _piece_length - _piece_length/2. - _detector_length/2. + _vertex_depth);

    vtx.setZ(Z);

    return vtx;
  }
  else {
    G4Exception("SciNearDetectorGeometry::GetVertex",
		"Geometry_S002", RunMustBeAborted,
		"ERROR.- Unknown region!");
  }
}



G4ThreeVector SciNearDetectorGeometry::RandomPointInBox(G4double width, 
							G4double height,
							G4double length,
							bool inVertex)
{
  G4double x=0, y=0, z=0;
  double m = (height/width - tan(atan(1.)/2.))/(1 - height/width*tan(atan(1.)/2.));
  if(IsOctagonal==1)
    do{
      x = (G4UniformRand() * width)  - width/2.;
      y = (G4UniformRand() * height) - height/2.;
      z = (G4UniformRand() * length) - length/2.;
    } while ( y >  width/2. * (m + tan(atan(1.)/2.)) + m*x  ||
	      y >  width/2. * (m + tan(atan(1.)/2.)) - m*x  ||
	      -y >  width/2. * (m + tan(atan(1.)/2.)) + m*x  ||
	      -y >  width/2. * (m + tan(atan(1.)/2.)) - m*x );
  else if(IsOctagonal==2)
    do{
      x = (G4UniformRand() * width)  - width/2.;
      y = (G4UniformRand() * height) - height/2.;
      z = -length/2.;
    } while ( (fabs(y) > height/2. * tan(atan(1.)/2)/2. + width/4. + x)  ||
	      (fabs(y) > height/2. * tan(atan(1.)/2)/2. + width/4. - x));
  else if(IsOctagonal==3)
    {
      x = (G4UniformRand() * width)  - width/2.;
      y = (G4UniformRand() * height) - height/2.;
      z = (G4UniformRand() * length) - length/2.;
    }
  else
    do{
      x = (G4UniformRand() * width)  - width/2.;
      y = (G4UniformRand() * height) - height/2.;
      z = (G4UniformRand() * length) - length/2.;
    } while (!inVertex && fabs(y) > height/2. * sqrt(1 - 4*x*x/pow(width,2))); 
  return G4ThreeVector(x, y, z);
}

G4ThreeVector SciNearDetectorGeometry::GaussianBeamSpot(G4ThreeVector origin,
						     G4double RMSx,
						     G4double RMSy){
  
  G4double x=0, y=0, z=0;
  G4double width = _piece_width, height = _piece_height;
  z = origin[2];
  double m = (height/width - tan(atan(1.)/2.))/(1 - height/width*tan(atan(1.)/2.));
  if(IsOctagonal==1){
    G4cout<<"Simulating beam spot in octagonal detector"<<G4endl; 
    do{
      x = G4RandGauss::shoot(origin[0], RMSx);
      y = G4RandGauss::shoot(origin[1], RMSy);
    } while ( y >  width/2. * (m + tan(atan(1.)/2.)) + m*x  ||
	      y >  width/2. * (m + tan(atan(1.)/2.)) - m*x  ||
	      -y >  width/2. * (m + tan(atan(1.)/2.)) + m*x  ||
	      -y >  width/2. * (m + tan(atan(1.)/2.)) - m*x );
  }
  else if(IsOctagonal==2){
    G4cout<<"Simulating beam spot in wide octagonal detector"<<G4endl; 
    do{
      x = G4RandGauss::shoot(origin[0], RMSx);
      y = G4RandGauss::shoot(origin[1], RMSy);
    } while ( (fabs(y) > height/2. * tan(atan(1.)/2)/2. + width/4. + x)  ||
	      (fabs(y) > height/2. * tan(atan(1.)/2)/2. + width/4. - x));
  }
  else
    do{
      x = G4RandGauss::shoot(origin[0], RMSx);
      y = G4RandGauss::shoot(origin[1], RMSy);
    } while ( fabs(y) > height/2. * sqrt(1 - 4*x*x/pow(width,2)));
  G4cout<<x<<"\t"<<y<<"\t"<<z<<G4endl;
  return G4ThreeVector(x,y,z);

}


void SciNearDetectorGeometry::CalculatePassiveTargetProb()
{
  const double Pi = 4*atan(1);
  G4int pass_A, act_A;  // number of nucleons (atomic number)
  G4double pass_M, act_M;  // molar mass
  G4double pass_dens, act_dens;  // density
  G4double pass_th = _passive_thickness;
  G4double act_th  = _number_active * _active_thickness;
  G4double pass_vol = _passive_thickness * 
    _piece_height * _piece_height * Pi /4.;
  G4double act_vol = _number_active * _active_thickness *
    _piece_height * _piece_height * Pi /4.;
  G4double vert_vol = _vertex_height * _vertex_width * _vertex_depth; 
  // tofix. some material properties can be fetched from G4Material
  act_A = 104;
  act_M = 104.149 * g/mole;
  act_dens = 1.06 * g/cm3;

  pass_A = 56;
  pass_M = 55.845 * g/mole;
  pass_dens = 7.874 * g/cm3;
  
  _passive_target_prob = (pass_A * pass_dens * pass_vol * act_M) /
    (pass_A * pass_dens * pass_vol * act_M + 
     act_A * act_dens * (act_vol + vert_vol) * pass_M);
}

G4int SciNearDetectorGeometry::GetRegion(G4double zpos)
{
  
  G4double npieces = (zpos + _detector_length/2)/_piece_length;
  
  G4double pieceFrac = npieces - (int)npieces;

  if ( pieceFrac < _passive_thickness ) return 1;

  return 0;
}






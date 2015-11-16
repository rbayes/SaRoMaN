// ----------------------------------------------------------------------------
///  \file   MindDetectorConstruction.h
///  \brief  User detector construction class.
/// 
///  \author   J Martin-Albo <jmalbos@ific.uv.es>
///  \date     14 Jun 2008
///  \version  $Id: MindDetectorConstruction.h 543 2014-11-01 23:03:10Z  $
///
///  Copyright (c) 2008, 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __DETECTOR_CONSTRUCTION__
#define __DETECTOR_CONSTRUCTION__

#include <G4VUserDetectorConstruction.hh>
#include <G4GDMLParser.hh>
#include "MindBarSD.h"
#include "MindFieldMapR.hh"
#include <G4SDManager.hh>
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4UserLimits.hh>


class G4VPhysicalVolume;
class SciNearDetectorGeometry;


/// TOFIX. Class description.
///
class MindDetectorConstruction: public G4VUserDetectorConstruction 
{
public:
  /// Constructor
  MindDetectorConstruction();
  /// Destructor
  ~MindDetectorConstruction();

  /// Returns the pointer to the physical volume that represents the world. 
  /// It is a Geant4-mandatory method, invoked by the G4RunManager before 
  /// starting the simulation run.
  G4VPhysicalVolume* Construct();
  // void ConstructSDandField();

  /// Returns a pointer to the detector geometry description
  SciNearDetectorGeometry* GetDetectorGeometry()
  { return _detector; }

private:
  void SetVolumeInformation(G4LogicalVolume* base, G4String detectorName);
  void SetAuxInformation(G4String basename, G4LogicalVolume* myvol,
			 const G4GDMLAuxListType auxlist);
  void SetNullField(G4LogicalVolume& detector_logic);
  void SetMagneticField(G4LogicalVolume& vol);
  SciNearDetectorGeometry* _detector;
  G4GDMLParser _gdml;
  std::string _gdml_file_name;
  bool _write_gdml;
  bool _use_gdml;
  std::vector<MindBarSD*> sensDetList;
  
};

#endif
  


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
  SciNearDetectorGeometry* _detector;
  std::string _gdml_file_name;
  bool _write_gdml;
  bool _use_gdml;
};

#endif
  


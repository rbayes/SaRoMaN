// ----------------------------------------------------------------------------
///  \file   SciNearDetectorGeometry.h
///  \brief  Definition of detector geometry.
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>
//   \date     14 Jun 2008
//   \version  $Id: SciNearDetectorGeometry.h 457 2011-09-07 17:20:21Z ryan $ 
//
//  Copyright (c) 2008, 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __DETECTOR_GEOMETRY__
#define __DETECTOR_GEOMETRY__

#include <G4ThreeVector.hh>

#include <vector>

class G4LogicalVolume;
class G4UniformMagField;


/// Definition of detector geometry.
//
class SciNearDetectorGeometry
{
public:
  /// Constructor
  SciNearDetectorGeometry();
  /// Destructor
  ~SciNearDetectorGeometry() {}
  
  /// Returns the logical volume representing the detector
  G4LogicalVolume* GetLogicalVolume() { return _logVolume; }

  /// Return probability of interaction in passive material
  G4double GetPassiveTargetProb() { return _passive_target_prob; }

  /// Returns a random vertex within 'region'
  G4ThreeVector GetVertex(G4String region);
  
  // Use a Gaussian beam spot with a pre-defined energy distribution
  G4ThreeVector GaussianBeamSpot(G4ThreeVector origin,
				 G4double RMSx,
				 G4double RMSy);

  /// Returns 'region' for a z position.
  G4int GetRegion(G4double zpos);
  
private:
  void Initialize();
  void SetInputParameters();
  G4LogicalVolume* DefineDetector();
  //void SetMagneticField();
  void SetNullField(G4LogicalVolume& detector_logic);
  void SetDipoleField(G4LogicalVolume& detector_logic);
  void SetMagneticField(G4LogicalVolume& detector_logic);
  G4ThreeVector RandomPointInBox(G4double, G4double, G4double, bool);
  void CalculatePassiveTargetProb();

private:

  //Set all volumes as sensitive?
  G4int _tasd;
  bool _isUniform;
  // Dimensions of PIECE and DETECTOR
  G4double _piece_width;
  G4double _piece_height;
  G4double _piece_length;
  G4double _mind_length;
  G4double _ear_width;
  G4double _ear_height;
  G4double _detector_length;
  G4double _bore_diameter;

  G4double _vertex_width;
  G4double _vertex_height;
  G4double _vertex_depth;
  bool _isVertexDet;
  G4double _TASDm_length;

  // Thicknesses of calorimeter layers
  G4double _active_thickness;
  G4double _passive_thickness;
  G4double _bracing_thickness;
  
  // Total number of PIECEs
  G4int _number_pieces;
  G4int _number_TASDpieces;

  // Number of scintillator planes per PIECE
  G4int _number_active;

  G4int IsOctagonal;

  // Gaps 1-3 are, respectively: separation between PASSIVE and 
  // 1st ACTIVE; separation between 1st and 2nd ACTIVE layers;
  // separation between 3rd layer and next PIECE.
  std::vector<G4double> _gaps;
  /* G4double _gap1; */
/*   G4double _gap2; */
/*   G4double _gap3; */

  // Displacement in Z-axis wrt PIECE system of reference (center)
  G4double _passive_displacement;
  std::vector<G4double> _active_displacement;
  std::vector<G4double> _brace_displacement;
  /* G4double _active1_displacement; */
/*   G4double _active2_displacement; */

  // Materials (name)
  G4String _active_material;
  G4String _passive_material;
  G4String _bracing_material;

  G4double _passive_target_prob;

  G4LogicalVolume* _logVolume;  ///< Pointer to detector logical volume

  G4UniformMagField* _magField; ///< Pointer to the magnetic field
};

#endif

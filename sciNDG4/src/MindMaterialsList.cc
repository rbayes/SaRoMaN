// ----------------------------------------------------------------------------
//  $Id: MindMaterialsList.cc 532 2013-05-21 14:47:07Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 4 Feb, 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindMaterialsList.h"

#include <G4NistManager.hh>
#include <G4Material.hh>



void MindMaterialsList::DefineMaterials()
{
  // Materials included in the Geant4 NIST database ...........................
  G4NistManager* nist = G4NistManager::Instance();
  
  // Air
  G4Material* G4_AIR = nist->FindOrBuildMaterial("G4_AIR");
  // Iron
  G4Material* G4_Fe = nist->FindOrBuildMaterial("G4_Fe");
  // Copper
  G4Material* G4_Cu = nist->FindOrBuildMaterial("G4_Cu");
  // Aluminum
  G4Material* G4_Al = nist->FindOrBuildMaterial("G4_Al");
  // Plastic Scintillator (polystyrene)
  G4Material* G4_POLYSTYRENE = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  // Vacuum (outer space)
  G4Material* G4_Galactic = nist->FindOrBuildMaterial("G4_Galactic");
}

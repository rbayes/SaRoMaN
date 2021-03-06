// ----------------------------------------------------------------------------
//  $Id: MindRunAction.cc 273 2009-06-09 15:29:07Z jmalbos $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindRunAction.h"
#include <G4Run.hh>



void MindRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
}



void MindRunAction::EndOfRunAction(const G4Run* run)
{
}


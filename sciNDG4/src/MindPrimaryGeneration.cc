// ----------------------------------------------------------------------------
//  $Id: MindPrimaryGeneration.cc 543 2014-11-01 23:03:10Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created:  Feb  8, 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindPrimaryGeneration.h"


MindPrimaryGeneration::MindPrimaryGeneration()
{
  Initialize();
}



void MindPrimaryGeneration::Initialize()
{
  G4String generator_name =
    MindConfigService::Instance().Generation().GetSParam("generator");
  if      (generator_name == "NUANCE") {
    _generator = new MindNuanceInterface(); 
  }
  else if (generator_name == "SINGLE_PARTICLE") {
    _generator = new MindSingleParticle();
  }
  else if (generator_name == "GENIE") {
    _generator = new MindGenieInterface();
  }
  else if (generator_name == "GENIE_EVENT") {
    _generator = new MindGenieEventInterface();
  }
  else if (generator_name == "CRY") {
    _generator = new MindCryInterface();
  }
  else {
  G4Exception("MindPrimaryGeneration::Initialize",
	      "MindPrimaryGeneration_S001", RunMustBeAborted,
	      "ERROR.- Unknown generator!");
  }
}



void MindPrimaryGeneration::GeneratePrimaries(G4Event* event)
{
  _generator->GeneratePrimaryVertex(event);
}




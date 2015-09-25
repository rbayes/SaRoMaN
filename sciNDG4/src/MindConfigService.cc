// ----------------------------------------------------------------------------
//  $Id: MindConfigService.cc 270 2009-06-04 14:44:06Z jmalbos $
//
//  Author :  J Martin-Albo <jmalbos@ific.uv.es>
//  Created:  16 Mar 2009
//      
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindConfigService.h"
#include "MindParamStore.h"

#include <bhep/sreader.h>
#include <bhep/gstore.h>



MindConfigService::MindConfigService():
  _geometry   (new MindParamStore("GEOMETRY")),
  _generation (new MindParamStore("GENERATION")),
  _physics    (new MindParamStore("PHYSICS")),
  _actions    (new MindParamStore("ACTIONS")),
  _job        (new MindParamStore("JOB"))
{
}



MindConfigService& MindConfigService::Instance()
{
  // The singleton instance is a static object:
  // it is initialized first time Instance() is invoked,
  // and deallocated when the program finishes.
  static MindConfigService svc;
  return svc;
}



void MindConfigService::SetConfigFile(const G4String& config_file)
{
  // Read general config file and fill corresponding stores
  _job       ->ParseFile(config_file);
  _geometry  ->ParseFile(config_file);
  _generation->ParseFile(config_file);
  _physics   ->ParseFile(config_file);
  _actions   ->ParseFile(config_file);
}

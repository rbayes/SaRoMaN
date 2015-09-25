// ----------------------------------------------------------------------------
//  $Id: MindParamStore.cc 543 2014-11-01 23:03:10Z  $
//  
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 29 May 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindParamStore.h"
#include <G4ExceptionHandler.hh>
#include <bhep/gstore.h>
#include <bhep/sreader.h>



MindParamStore::MindParamStore(const G4String& keyword):
  _store(new bhep::gstore()), _keyword(keyword)
{
}
  
  

MindParamStore::~MindParamStore()
{
  delete _store;
}
  
  
  
void MindParamStore::ParseFile(const G4String& filename)
{
  bhep::sreader parser(*_store);
  parser.file(filename);
  parser.info_level(bhep::MUTE);
  parser.group(_keyword);
  parser.read();
}
  
  
  
G4int MindParamStore::GetIParam(const G4String& name) const
{
  try {
    return _store->fetch_istore(name);
  }
  catch (...) {
    G4Exception("MindParamStore::GetIParam",
		"MindParamStore_S001", RunMustBeAborted,
		"ERROR.- Fetch failed!");
  }
}

G4double MindParamStore::GetDParam(const G4String& name) const
{
  try {
    return _store->fetch_dstore(name);
  }
  catch (...) {
    G4Exception("MindParamStore::GetDParam",
		"MindParamStore_S002", RunMustBeAborted,
		"ERROR.- Fetch failed!");
  }
}

G4String MindParamStore::GetSParam(const G4String& name) const
{
  try {
    return _store->fetch_sstore(name);
  }
  catch (...) {
    G4Exception("MindParamStore::GetSParam",
		"MindParamStore_S003", RunMustBeAborted,
		"ERROR.- Fetch failed!");
  }
}
  
const std::vector<G4double> MindParamStore::GetVParam(const G4String& name) const
{
  try {
    return _store->fetch_vstore(name);
  }
  catch (...) {
    G4Exception("MindParamStore::GetVParam",
		"MindParamStore_S004", RunMustBeAborted,
		"ERROR.- Fetch failed!");
  }
}



G4bool MindParamStore::PeekIParam(const G4String& name) const
{
  return _store->find_istore(name);
}

G4bool MindParamStore::PeekDParam(const G4String& name) const
{
  return _store->find_dstore(name);
}

G4bool MindParamStore::PeekSParam(const G4String& name) const
{
  return _store->find_sstore(name);
}

G4bool MindParamStore::PeekVParam(const G4String& name) const
{
  return _store->find_vstore(name);
}



void MindParamStore::Error(const G4String& name) const
{
  G4Exception("MindParamStore::Error",
	      "MindParamStore_S005", RunMustBeAborted,
	      "ERROR.- Param not found in store.");
}

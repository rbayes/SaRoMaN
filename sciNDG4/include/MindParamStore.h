// ----------------------------------------------------------------------------
///  \file   MindParamStore.h
///  \brief  A collection of configuration parameters.
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>    
///  \date     29 May 2009
///  \version  $Id: MindParamStore.h 543 2014-11-01 23:03:10Z  $
///
///  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __PARAM_STORE__
#define __PARAM_STORE__

#include <globals.hh>
#include <vector>
#include <CLHEP/Units/GlobalSystemOfUnits.h>

namespace bhep { class gstore; }


/// A collection of configuration parameters.
///
class MindParamStore
{
public:
  /// Constructor
  MindParamStore(const G4String& groupID);
  /// Destructor
  ~MindParamStore();
  
  /// Return an integer from the parameter store
  G4int GetIParam(const G4String& name) const;
  /// Return a double from the parameter store
  G4double GetDParam(const G4String& name) const;
  /// Return a string from the parameter
  G4String GetSParam(const G4String& name) const;
  /// Return a vdouble from the store
  const std::vector<G4double> GetVParam(const G4String&) const;
  
  /// Parse file and fill store
  void ParseFile(const G4String&);
  
  // Check whether a parameters exists or not in the store
  //
  /// Check integer-parameter existence in store
  G4bool PeekIParam(const G4String& name) const;
  /// Check double-parameter existence in store
  G4bool PeekDParam(const G4String& name) const;
  /// Check string-parameter existence in store
  G4bool PeekSParam(const G4String& name) const;
  /// Check vector double-parameter existence in store
  G4bool PeekVParam(const G4String& name) const;


private:
  void Error(const G4String& name) const;
  
private:
  bhep::gstore* _store;
  G4String _keyword;
};

#endif

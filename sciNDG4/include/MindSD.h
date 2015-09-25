// ----------------------------------------------------------------------------
///  \file   MindSD.h
///  \brief  
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>
///  \date     15 Apr 2009
///  \version  $Id: MindSD.h 280 2009-06-11 12:36:45Z jmalbos $
///
///  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#ifndef __SENSITIVE_DETECTOR__
#define __SENSITIVE_DETECTOR__

#include <G4VSensitiveDetector.hh>
#include "MindHit.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


/// TOFIX. 
///

class MindSD: public G4VSensitiveDetector
{
public:
  /// Constructor
  MindSD(G4String);
  /// Destructor
  ~MindSD() {}

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*) {}

private:
  MindHitsCollection* _MHCollection;
};

#endif

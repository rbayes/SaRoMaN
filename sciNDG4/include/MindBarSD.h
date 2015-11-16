// ----------------------------------------------------------------------------
///  \file   MindBarSD.h
///  \brief  
///
///  \author   R Bayes (ryan.bary@glasgow.ac.uk
///  \date     15 Apr 2009
///  \version  $Id: MindBarSD.h 280 2009-06-11 12:36:45Z ryan $
///
///  Copyright (c) 2015 Glasgow Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __SENSITIVE_DETECTOR__
#define __SENSITIVE_DETECTOR__

#include <G4VSensitiveDetector.hh>
#include "MindBarHit.h"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;



class MindBarSD: public G4VSensitiveDetector
{
public:
  /// Constructor
  MindBarSD(G4String);
  /// Destructor
  ~MindBarSD() {}

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*) {}

private:
  MindBarHitsCollection* _MHCollection;
  std::vector<G4String> _volumelist;
};

#endif

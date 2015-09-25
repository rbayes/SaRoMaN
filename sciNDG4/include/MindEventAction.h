// ----------------------------------------------------------------------------
///  \file   MindEventAction.h
///  \brief  User action at event scope.
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>
///  \date     15 Apr 2009
///  \version  $Id: MindEventAction.h 462 2011-09-20 15:41:28Z ryan $ 
///
///  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#ifndef __EVENT_ACTION__
#define __EVENT_ACTION__

#include <G4UserEventAction.hh>

#include <bhep/hit.h>

class G4Event;
class G4HCofThisEvent;
//class EventManager2;
//class writer_root;


/// TOFIX. Class description.
///

class MindEventAction: public G4UserEventAction
{
public:
  /// Constructor
  MindEventAction();
  /// Destructor
  ~MindEventAction();
  
  ///
  void BeginOfEventAction(const G4Event*);
  
  ///
  void EndOfEventAction(const G4Event*);

private:
  void Initialize();
  void ProcessHits(G4HCofThisEvent*);

  bhep::particle* _leptonShowerPart;
  bhep::particle* _hadronShowerPart;

  int _evtNo;
};

#endif

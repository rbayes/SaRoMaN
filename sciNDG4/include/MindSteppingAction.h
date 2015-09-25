// ----------------------------------------------------------------------------
///  \file   MindSteppingAction.h
///  \brief
///
///  \author   A Laing <a.laing@physics.gla.ac.uk>
///  \date     3 Jun 2009
///  \version  $Id: MindSteppingAction.h 331 2009-11-05 10:08:16Z alaing $
///
///  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __STEPPING_ACTION__
#define __STEPPING_ACTION__

#include <G4UserSteppingAction.hh>
#include <G4ThreeVector.hh>

class G4Step;

class MindSteppingAction: public G4UserSteppingAction
{
 public:
  /// Constructor
  MindSteppingAction()  {}
  /// Destructor
  ~MindSteppingAction() {}


  /// The action
  void UserSteppingAction(const G4Step*);

  ///So the check can be reset for each track.
  void reset_stepCheck(){ _vertCheck = false; }

 private:
  ///used to check very long tracks which might be looping at endpoint.
  bool _vertCheck;
  G4ThreeVector _lastPoint;//< position at last check.

};

#endif

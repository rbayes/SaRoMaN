// ----------------------------------------------------------------------------
///  \file   MindTrackingAction.h
///  \brief
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>
///  \date     3 Jun 2009
///  \version  $Id: MindTrackingAction.h 352 2009-09-22 13:36:52Z alaing $
///
///  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __TRACKING_ACTION__
#define __TRACKING_ACTION__

#include <G4UserTrackingAction.hh>
#include <G4VProcess.hh>

class G4Track;


class MindTrackingAction: public G4UserTrackingAction
{
 public:
  /// Constructor
  MindTrackingAction()  {}
  /// Destructor
  ~MindTrackingAction() {}
  
  /// Operations before the processing of a track
  void PreUserTrackingAction(const G4Track*);
  /// Operation after the processing of a track
  void PostUserTrackingAction(const G4Track*);
  
 private:
  //Stores track id of primary lepton.
  G4int _primLep;

};

#endif

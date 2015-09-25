// ----------------------------------------------------------------------------
//  $Id: MindSteppingAction.cc 371 2009-11-03 09:30:35Z alaing $
//
//  Author : A Laing <a.laing@physics.gla.ac.uk>
//  Created: 5 October 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindSteppingAction.h"

#include <G4Step.hh>
#include <G4Track.hh>

void MindSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();

  // To protect against inifinite loops.
  // Kill any particle with more than 1000 steps that is within
  // 1 mm displacement from its origin or
  // one that moves less than 1mm in 500 steps thereafter.
  G4int stepNo = track->GetCurrentStepNumber();
  
  if ( stepNo > 1000 ) {
    
    if ( !_vertCheck ){
      G4ThreeVector disp = track->GetPosition()
	- track->GetVertexPosition();
      
      if ( disp.mag() <= 1.0 ){
	G4cout << "GOING TO KILL INFINITE LOOP TRACK" << G4endl;
	track->SetTrackStatus( fStopAndKill );
      } else {
	_vertCheck = true;
	_lastPoint = track->GetPosition();
      }
    } else if ( stepNo % 500 == 0 ){
      G4ThreeVector disp = track->GetPosition()
	- _lastPoint;
      
      if ( disp.mag() <= 1.0 ){
	G4cout << "Killing track looping at endpoint" << G4endl;
	track->SetTrackStatus( fStopAndKill );
      } else
	_lastPoint = track->GetPosition();
      
    }
    
  }
  
}

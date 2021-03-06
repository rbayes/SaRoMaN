// ----------------------------------------------------------------------------
//  $Id: MindHit.cc 532 2013-05-21 14:47:07Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindHit.h"
#include <G4Allocator.hh>


G4Allocator<MindHit> MindHitAllocator;



MindHit::MindHit(const MindHit& right): G4VHit()
{
  _track_id   = right._track_id;
  _position   = right._position;
  _energy_dep = right._energy_dep;
  _time       = right._time;
}
  


const MindHit& MindHit::operator=(const MindHit& right)
{
  _track_id   = right._track_id;
  _position   = right._position;
  _energy_dep = right._energy_dep;
  _time       = right._time;
  return *this;
}



G4int MindHit::operator==(const MindHit& right) const
{
  return (this==&right) ? 1 : 0;
}

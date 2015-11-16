// ----------------------------------------------------------------------------
//  $Id: MindHit.cc 532 2013-05-21 14:47:07Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindBarHit.h"
#include <G4Allocator.hh>


G4Allocator<MindBarHit> MindBarHitAllocator;



MindBarHit::MindBarHit(const MindBarHit& right): G4VHit()
{
  _track_id   = right._track_id;
  _position   = right._position;
  _energy_dep = right._energy_dep;
  _time       = right._time;
  _barCopy    = right._barCopy;
  _module     = right._module;
  _isYbar     = right._isYbar;
}
  


const MindBarHit& MindBarHit::operator=(const MindBarHit& right)
{
  _track_id   = right._track_id;
  _position   = right._position;
  _energy_dep = right._energy_dep;
  _time       = right._time;
  _barCopy    = right._barCopy;
  _module     = right._module;
  _isYbar     = right._isYbar;
  return *this;
}



G4int MindBarHit::operator==(const MindBarHit& right) const
{
  return (this==&right) ? 1 : 0;
}

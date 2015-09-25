// ----------------------------------------------------------------------------
//  $Id: MindHit.h 532 2013-05-21 14:47:07Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#ifndef __HIT__
#define __HIT__

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>


class MindHit: public G4VHit
{
public:
  /// Constructor
  MindHit() {}
  /// Destructor
  ~MindHit() {}
  /// Copy-constructor
  MindHit(const MindHit&);
  /// Assignement operator
  const MindHit& operator=(const MindHit&);
  /// Equality operator
  G4int operator==(const MindHit&) const;

  /// Memory allocation
  inline void* operator new(size_t);
  inline void  operator delete(void*);

public:
  /// Setters & Getters
  ///
  G4int GetTrackID() { return _track_id; }
  void  SetTrackID(G4int id) { _track_id = id; }
  ///
  G4ThreeVector GetPosition() { return _position; }
  void SetPosition(G4ThreeVector p) { _position = p; }
  ///
  G4double GetEnergyDeposit() { return _energy_dep; }
  void SetEnergyDeposit(G4double e) { _energy_dep = e; }
  //
  G4double GetHitTime() { return _time; }
  void SetHitTime(G4double t) { _time = t; }

private:
  G4int _track_id;
  G4double _energy_dep;
  G4ThreeVector _position;
  G4double _time;
};


typedef G4THitsCollection<MindHit> MindHitsCollection;


extern G4Allocator<MindHit> MindHitAllocator;


inline void* MindHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) MindHitAllocator.MallocSingle();
  return aHit;
}


inline void MindHit::operator delete(void *aHit)
{
  MindHitAllocator.FreeSingle((MindHit*) aHit);
}


#endif

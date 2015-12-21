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


class MindBarHit: public G4VHit
{
public:
  /// Constructor
  MindBarHit() {}
  /// Destructor
  ~MindBarHit() {}
  /// Copy-constructor
  MindBarHit(const MindBarHit&);
  /// Assignement operator
  const MindBarHit& operator=(const MindBarHit&);
  /// Equality operator
  G4int operator==(const MindBarHit&) const;

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

  G4ThreeVector GetBarTranslation() { return _barCopy; }
  void SetBarTranslation(G4ThreeVector b) { _barCopy = b; }

  G4String GetModule() { return _module; }
  void SetModule(G4String m) { _module = m; }

  G4int GetBarOrientation() { return _isYbar; }
  void SetBarOrientation(G4int i) { _isYbar = i; }

  G4int GetBarNumber() { return _barNumber; } 
  void SetBarNumber(G4int i) { _barNumber = i; }
private:
  G4int _track_id;
  G4double _energy_dep;
  G4ThreeVector _position;
  G4ThreeVector _barCopy;
  G4double _time;
  G4String _module;
  G4int _isYbar;
  G4int _barNumber;
};


typedef G4THitsCollection<MindBarHit> MindBarHitsCollection;


extern G4Allocator<MindBarHit> MindBarHitAllocator;


inline void* MindBarHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) MindBarHitAllocator.MallocSingle();
  return aHit;
}


inline void MindBarHit::operator delete(void *aHit)
{
  MindBarHitAllocator.FreeSingle((MindBarHit*) aHit);
}


#endif

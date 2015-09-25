#include "MindLookupTable.h"

MindLookupTable::MindLookupTable()
{
  clear();
}

MindLookupTable& MindLookupTable::Instance()
{
  // The singleton instance is a static object:
  // it is initialized first time Instance() is invoked,
  // and deallocated when the program finishes.
  static MindLookupTable svc;
  return svc;
}

void MindLookupTable::clear()
{
  //Clear vectors.
  _primaries.clear();

  _leptonShower.clear();

  _hadronShower.clear();

  _primlep = false;

}

void MindLookupTable::add_particle(G4int shower, G4int particle)
{
  //Add a particle identity to the appropriate shower vector.

  if ( shower == 1 )
    _leptonShower[particle] = true;

  else if ( shower == 2 )
    _hadronShower[particle] = true;

  else
    G4Exception("LookupTable::add_particle",
		"Lookup_S001", RunMustBeAborted,
		"ERROR.- Incorrect shower number!");

}

void MindLookupTable::add_primary(G4int particle)
{
  _primaries[particle] = true;
}

G4int MindLookupTable::find_particle(G4int particle)
{
  //Check if particle appears in one of the showers
  //return 0 if not present, 1 for lepton and 2 for hadron.

  std::map<G4int,bool>::iterator it1;
  
  it1 = _leptonShower.find( particle );
  if ( it1 != _leptonShower.end() ) return 1;

  it1 = _hadronShower.find( particle );
  if ( it1 != _hadronShower.end() ) return 2;

  return 0;
}

bool MindLookupTable::lepton_shower()
{
  //return true if particles appended to lepton vector.

  if ( _leptonShower.size() == 0 ) return false;

  return true;
}

bool MindLookupTable::hadron_shower()
{
  //return true if particles appended to hadron vector.

  if ( _hadronShower.size() == 0 ) return false;

  return true;
}

bool MindLookupTable::check_primary(G4int particle)
{
  //return true if id is in the primaries vector.

  std::map<G4int,bool>::iterator it2;

  it2 = _primaries.find( particle );
  if ( it2 != _primaries.end() ) return true;
  
  return false;
}

void MindLookupTable::set_primLep(G4int PDG, G4double Px, G4double Py, G4double Pz)
{

  _primlep = true;
  
  _primLepPDG = PDG;

  _primLepStartMom = G4ThreeVector( Px, Py, Pz );

}

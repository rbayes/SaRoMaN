//Static class to contain the identifiers of those particles to be put in
//the shower particles.
//Can be extended to store all particle details for output at the end of the
//event.
//Andrew Laing, 08/09.

#ifndef __MINDLOOKUPTABLE__
#define __MINDLOOKUPTABLE__

#include <map>
#include <globals.hh>

#include <G4ThreeVector.hh>

class MindLookupTable
{
 public:
  /// Returns the singleton class instance
  static MindLookupTable& Instance();
  //Add a particle, shower=1 lepton, =2 hadron.
  void add_particle(G4int shower, G4int particle);
  //Add a primary id to the primaries vector.
  void add_primary(G4int particle);
  //return 0 if particle not found, 1 if in lepton shower, 2 in hadron.
  G4int find_particle(G4int particle);
  //Clear out vectors.
  void clear();
  //Check if lepton shower has particles.
  bool lepton_shower();
  //check if hadron shower has particles.
  bool hadron_shower();
  //Check if particle id is a primary.
  bool check_primary(G4int particle);
  //Check if a primary lepton is set.
  bool check_primLep(){ return _primlep; }
  //set primary lepton.
  void set_primLep(G4int PDG, G4double Px, G4double Py, G4double Pz);
  //Get primary lepton info.
  G4int get_primLepPDG(){ return _primLepPDG; }
  G4ThreeVector get_pLMom(){ return _primLepStartMom; }

 private:
  /// constructor (hidden)
  MindLookupTable();
  /// destructor (hidden)
  ~MindLookupTable() {}
  /// copy-constructor (hidden)
  MindLookupTable(MindLookupTable const&);
  /// assignement operator (hidden)
  MindLookupTable& operator=(MindLookupTable const&);

 private:

  std::map<G4int,bool> _primaries;
  std::map<G4int,bool> _leptonShower;
  std::map<G4int,bool> _hadronShower;

  bool _primlep; //True if a primary lepton is being tracked (i.e not decayed by genie)
  G4int _primLepPDG; //Set if _primLep true.
  G4ThreeVector _primLepStartMom; //Momentum from Generator (should be same at PreTrack)
  
};

#endif

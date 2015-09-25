// ----------------------------------------------------------------------------
//  $Id: MindSD.cc 543 2014-11-01 23:03:10Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindSD.h"
#include "MindHit.h"

#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>



MindSD::MindSD(G4String name): G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="MindCollection");
}



void MindSD::Initialize(G4HCofThisEvent* HCE)
{
  _MHCollection = 
    new MindHitsCollection(SensitiveDetectorName, collectionName[0]);

  static G4int HCID = -1;
  if (HCID < 0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  
  HCE->AddHitsCollection(HCID, _MHCollection);
}



G4bool MindSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  
  G4double edep = step->GetTotalEnergyDeposit();
  G4double time = step->GetDeltaTime();
  if (edep == 0.) return false;
  MindHit* hit = new MindHit();
  hit->SetTrackID(step->GetTrack()->GetTrackID());
  hit->SetPosition(step->GetPostStepPoint()->GetPosition());
  // G4cout<<edep<<"\t"<<step->GetPostStepPoint()->GetPosition()[2]<<G4endl;
  hit->SetEnergyDeposit(edep);
  hit->SetHitTime(time);
  _MHCollection->insert(hit);
}



// void MindSD::EndOfEvent(G4HCofThisEvent*)
// {

// }



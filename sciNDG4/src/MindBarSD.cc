// ----------------------------------------------------------------------------
//  $Id: MindSD.cc 543 2014-11-01 23:03:10Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindBarSD.h"
#include "MindBarHit.h"

#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>



MindBarSD::MindBarSD(G4String name): G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="MindCollection");
  
  std::string delimiter = "/";
  
  size_t pos = 0;
  std::string token;
  std::string s = name;
  while ((pos = s.find(delimiter)) != std::string::npos) { 
    token = s.substr(0, pos);
    _volumelist.push_back(token);
    s.erase(0, pos + delimiter.length());
    // std::cout<<token<<std::endl;
  }
  // std::cout<<s<<std::endl;
  _volumelist.push_back(s);
}



void MindBarSD::Initialize(G4HCofThisEvent* HCE)
{
  _MHCollection = 
    new MindBarHitsCollection(SensitiveDetectorName, collectionName[0]);

  static G4int HCID = -1;
  if (HCID < 0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  
  HCE->AddHitsCollection(HCID, _MHCollection);
}



G4bool MindBarSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();

  // G4Solid* modulesolid = theTouchable->GetVolume(1);
  // Get the z-dimension of the solid
  
  // if (_volumelist[2].contains
  // Record plane z position as part of the bar translation
    

  G4double edep = step->GetTotalEnergyDeposit();
  G4double time = step->GetDeltaTime();
  if (edep == 0.) return false;
  MindBarHit* hit = new MindBarHit();
  hit->SetTrackID(step->GetTrack()->GetTrackID());
  hit->SetPosition(step->GetPostStepPoint()->GetPosition());
  // G4cout<<edep<<"\t"<<step->GetPostStepPoint()->GetPosition()[2]<<G4endl;
  hit->SetEnergyDeposit(edep);
  hit->SetHitTime(time);

  G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector copyTrans = theTouchable->GetTranslation();
  /*
  G4ThreeVector barTrans = theTouchable->GetTranslation(1);
  G4ThreeVector modTrans = theTouchable->GetTranslation(2);
  G4ThreeVector detTrans = theTouchable->GetTranslation(3);
  std::cout<<" hit position = ("<<pos[0]
	   <<","<<pos[1]<<","<<pos[2]<<")\n";
  std::cout<<" Copy Translation = ("<<copyTrans[0]
	   <<","<<barTrans[1]<<","<<copyTrans[2]<<")\n";
  std::cout<<" Bar Translation = ("<<barTrans[0]
	   <<","<<barTrans[1]<<","<<barTrans[2]<<")\n\n";
  std::cout<<" Module Translation = ("<<modTrans[0]
	   <<","<<modTrans[1]<<","<<modTrans[2]<<")\n";
  std::cout<<" Detector Translation = ("<<detTrans[0]
	   <<","<<detTrans[1]<<","<<detTrans[2]<<")\n\n";
  */
  G4double barOffset = 0;

  hit->SetBarTranslation(copyTrans);
  hit->SetModule(_volumelist[2]);
  if(_volumelist.back().contains("Y")){
    hit->SetBarOrientation(1);
  } else {
    hit->SetBarOrientation(0);
  }

  _MHCollection->insert(hit);
}

// void MindSD::EndOfEvent(G4HCofThisEvent*)
// {

// }


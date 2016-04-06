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
  
  std::cout<<collectionName.size()<<std::endl;
  
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
  //G4double time = step->GetDeltaTime();
  //G4double time = step->GetPostStepPoint()->GetGlobalTime();
  //G4double time = step->GetPostStepPoint()->GetLocalTime();
  //G4double time = step->GetPostStepPoint()->GetProperTime();
  //G4double time = step->GetTrack()->GetProperTime(); 
  G4double time = step->GetTrack()->GetGlobalTime();
  if (edep == 0.) return false;
  MindBarHit* hit = new MindBarHit();
  hit->SetTrackID(step->GetTrack()->GetTrackID());
  hit->SetPosition(step->GetPostStepPoint()->GetPosition());
  // G4cout<<edep<<"\t"<<step->GetPostStepPoint()->GetPosition()[2]<<G4endl;
  hit->SetEnergyDeposit(edep);
  hit->SetHitTime(time);

  G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector copyTrans = theTouchable->GetTranslation();
  G4String volName = theTouchable->GetVolume()->GetName();
  G4int copyNo = theTouchable->GetCopyNumber();
  
  
  std::cout << volName;
  std::cout<<" hit position = ("<<pos.x()
  	   <<","<<pos.y()<<","<<pos.z()<<"), copyNo = "<<copyNo;
  std::cout<<" Copy Translation = ("<<copyTrans.x()
  	   <<","<<copyTrans.y()<<","<<copyTrans.z()<<")\n";
  
  
  G4double barOffset = 0;

  hit->SetBarTranslation(copyTrans);
  
  // hit->SetModule(_volumelist[2]);
  if(volName.contains("BarY")){
    // std::cout << volName << std::endl;
    hit->SetBarOrientation(1);
  } else {
    hit->SetBarOrientation(0);
  }
  hit->SetBarNumber(copyNo);

  _MHCollection->insert(hit);
}

// void MindSD::EndOfEvent(G4HCofThisEvent*)
// {

// }



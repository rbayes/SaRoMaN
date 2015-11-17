// ----------------------------------------------------------------------------
//  $Id: MindEventAction.cc 543 2014-11-01 23:03:10Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 15 Apr 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group 
// ----------------------------------------------------------------------------

#include "MindEventAction.h"
#include "MindConfigService.h"
#include "MindParamStore.h"
#include "MindBarSD.h"
#include "MindUtils.h"
#include "MindLookupTable.h"

#include <G4Event.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4VVisManager.hh>
#include <G4SDManager.hh>
#include <G4UImanager.hh> 

#include <bhep/bhep_svc.h>
//#include <bhep/EventManager2.h>
//#include <bhep/writer_root.h>

using namespace bhep;

MindEventAction::MindEventAction():
  _evtNo(0)
{
}


MindEventAction::~MindEventAction()
{
}



void MindEventAction::BeginOfEventAction(const G4Event* event)
{
  // retrieve bhep transient event and clear
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();
  
  // set Geant4 event id property
  bevt.set_event_number(_evtNo);
  G4int eventID = event->GetEventID();
  
  G4cout << "Event: " << eventID << G4endl;
  
  bevt.add_property("G4EventID", eventID);
}



void MindEventAction::EndOfEventAction(const G4Event* event)
{
  // Visualization of tracks ........................................
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  
  G4int number_trajectories;
  if (trajectoryContainer)
    number_trajectories = trajectoryContainer->entries();
  
  if (G4VVisManager::GetConcreteInstance()) {
    for (G4int i=0; i<number_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory*) ((*(event->GetTrajectoryContainer()))[i]);
      trj->DrawTrajectory();
    }
  }
  // Make bhep containers for the hits associated to particles without bhep particles.
  if ( MindLookupTable::Instance().lepton_shower() ){
    _leptonShowerPart = new bhep::particle( bhep::TRUTH, "lepton_shower");
    _leptonShowerPart->add_property("length", 0.0);
    _leptonShowerPart->add_property("CreatorProcess", "showering"); }

  if ( MindLookupTable::Instance().hadron_shower() ){
    _hadronShowerPart = new bhep::particle( bhep::TRUTH, "hadron_shower");
    _hadronShowerPart->add_property("length", 0.0);
    _hadronShowerPart->add_property("CreatorProcess", "showering"); }

  // Hits in sensitive detectors ....................................
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  ProcessHits(HCE);
  
  // bhep event is written in dst
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();
  
  if ( MindLookupTable::Instance().lepton_shower() ) bevt.add_true_particle( _leptonShowerPart );
  if ( MindLookupTable::Instance().hadron_shower() ) bevt.add_true_particle( _hadronShowerPart );
  //
  
  bhep::bhep_svc::instance()->get_writer_root().write(bevt, _evtNo);
  _evtNo++;
  bevt.clear();
  
  MindLookupTable::Instance().clear();
  
}



void MindEventAction::ProcessHits(G4HCofThisEvent* HCE)
{
  G4int ID = -1; //dummy value to start.
  G4int pstatus;

  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();
  bhep::particle* bpart;
  
  G4SDManager* SDman  = G4SDManager::GetSDMpointer();
  G4int collection_id = SDman->GetCollectionID("MindCollection");
  
  MindHitsCollection* THC =
    (MindHitsCollection*)(HCE->GetHC(collection_id));
  
  for (G4int i=0; i<(THC->entries()); i++) {
  
    bhep::hit* bhit = new bhep::hit("tracking");
  
    G4ThreeVector xyz = (*THC)[i]->GetPosition();
    // G4cout<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\n";
    bhit->set_point(bhep::Point3D(xyz.x(), xyz.y(), xyz.z()));
    
    G4double energy_dep = (*THC)[i]->GetEnergyDeposit();
    bhit->add_property("EnergyDep", energy_dep);
    
    G4double time = (*THC)[i]->GetHitTime();
    bhit->add_property("time", time);

    G4double transbarpos = (*THC)[i]->GetBarOrientation()==1 ?
      (*THC)[i]->GetBarTranslation()[1] : (*THC)[i]->GetBarTranslation()[0];
    bhit->add_property("barPosT", transbarpos);

    G4double longbarpos = (*THC)[i]->GetBarTranslation()[2];
    bhit->add_property("barPosZ", longbarpos);

    G4int barorientation = (*THC)[i]->GetBarOrientation();
    bhit->add_property("IsYBar", barorientation);

    G4String module = (*THC)[i]->GetModule();
    bhit->add_property("detmodule", module);
    
    pstatus = MindLookupTable::Instance().find_particle( (*THC)[i]->GetTrackID() );

    if ( pstatus == 0 ){
    
      if ( ID != (*THC)[i]->GetTrackID() ){
	ID = (*THC)[i]->GetTrackID();
	bpart = MindUtils::GetBParticle(ID);
      }

      bhit->set_mother_particle(*bpart);
      bpart->add_hit("tracking", bhit);

    } else if ( pstatus == 1 ){

      bhit->set_mother_particle(*_leptonShowerPart);
      _leptonShowerPart->add_hit("tracking", bhit);

    } else {

      bhit->set_mother_particle(*_hadronShowerPart);
      _hadronShowerPart->add_hit("tracking", bhit);

    }

  }
}

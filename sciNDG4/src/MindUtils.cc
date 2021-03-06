// ----------------------------------------------------------------------------
//  $Id: MindUtils.cc 543 2014-11-01 23:03:10Z  $
//  
//  Author :  J Martin-Albo <jmalbos@ific.uv.es>
//  Created:  10 June 2009
//
//  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindUtils.h"
#include <bhep/bhep_svc.h>
#include <bhep/particle.h>


bhep::particle* MindUtils::GetBParticle(G4int G4TrackID)
{
  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();
  vector<bhep::particle*> bparts;
  
  if (bevt.filter(bhep::TRUTH, "G4TrackID", bhep::to_string(G4TrackID), bparts))
    return bparts[0];
  else {
    G4Exception("MindUtils::GetBParticle",
		"Utils_S001", RunMustBeAborted,
		"ERROR.- Particle does not exist in bhep transient event!");
  }
}

void MindUtils::OpenOutputDst(const G4String& dst_name, G4int seed_val)
{
  if (!bhep::bhep_svc::instance()->get_writer_root().isOpen()) {
    
    bhep::bhep_svc::instance()->get_writer_root().open(dst_name, "RECREATE");
    
    // set dst info
    time_t time1; time(&time1);
    bhep::bhep_svc::instance()->get_dst_info().set_time(time1);
    bhep::bhep_svc::instance()->get_dst_info().set_type(bhep::MC);
    bhep::bhep_svc::instance()->get_dst_info().set_label("G4outputDST");
    
    // set run info
    bhep::bhep_svc::instance()->get_run_info().set_time(time1);
    bhep::bhep_svc::instance()->get_run_info().set_type(bhep::MC);
    bhep::bhep_svc::instance()->get_run_info().set_label("G4run");

    if ( seed_val < 0 )
      bhep::bhep_svc::instance()->get_run_info().add_property( "seed", (G4int)time(0) );
    else
      bhep::bhep_svc::instance()->get_run_info().add_property( "seed", seed_val );
  }
}

void MindUtils::CloseOutputDst()
{ 
  bhep::dst& ndst = bhep::bhep_svc::instance()->get_dst_info();
  bhep::run& nrun = bhep::bhep_svc::instance()->get_run_info();
  
  // add final parameters to dst info
  size_t nevt = bhep::bhep_svc::instance()->get_writer_root().max_events();
  ndst.set_Nevt(nevt);
  
  // save dst info
  bhep::bhep_svc::instance()->get_writer_root().write_dst_info(ndst);
  bhep::bhep_svc::instance()->get_dst_info().clear();
  
  // save run info
  bhep::bhep_svc::instance()->get_writer_root().write_run_info(nrun);
  bhep::bhep_svc::instance()->get_run_info().clear();
  
  // close file
  bhep::bhep_svc::instance()->get_writer_root().close();
}


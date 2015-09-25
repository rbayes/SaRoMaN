#include "WriteUtil.h"
#include <bhep/bhep_svc.h>


void WriteUtil::OpenOutputDst(const std::string& dst_name)
{
  if (!bhep::bhep_svc::instance()->get_writer_root().isOpen()) {
    
    bhep::bhep_svc::instance()->get_writer_root().open(dst_name, "RECREATE");

    // set dst info
    time_t time1; time(&time1);
    bhep::bhep_svc::instance()->get_dst_info().set_time(time1);
    bhep::bhep_svc::instance()->get_dst_info().set_type(bhep::MC);
    bhep::bhep_svc::instance()->get_dst_info().set_label("DigioutputDST");
    
    // set run info
    bhep::bhep_svc::instance()->get_run_info().set_time(time1);
    bhep::bhep_svc::instance()->get_run_info().set_type(bhep::MC);
    bhep::bhep_svc::instance()->get_run_info().set_label("Digirun");
  }
}


void WriteUtil::CloseOutputDst()
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

#include <MINDfitman.h>
#include <mind/MINDsetup.h>
#include <mind/SetupSk.h>

MINDfitman::MINDfitman()
{
}

MINDfitman& MINDfitman::instance()
{

  static MINDfitman svc;
  return svc;

}

void MINDfitman::set_man_parameters(const bhep::gstore params, Setup& det)
{
  //set the parameters which are common to both modes
  //and attach the parameter store.
  _store = params;

  config_common_props();

  config_geometry( det );

  config_fitter();

  config_navigator();

  config_model();

  config_matching();

}

void MINDfitman::fit_mode()
{

  _man.fitting_svc().retrieve_fitter<KalmanFitter>(_fitter,_model)
    .set_max_local_chi2ndf( _store.fetch_dstore("chi2node_max") );

  _man.fitting_svc().retrieve_fitter<KalmanFitter>(_fitter,_model)
    .set_number_allowed_outliers( _store.fetch_istore("max_outliers") );
  
  _inRecMode = false;

}

void MINDfitman::rec_mode()
{

  _man.fitting_svc().retrieve_fitter<KalmanFitter>(_fitter,_model)
    .set_max_local_chi2ndf( _store.fetch_dstore("pat_rec_max_chi") );

  _man.fitting_svc().retrieve_fitter<KalmanFitter>(_fitter,_model)
    .set_number_allowed_outliers( _store.fetch_istore("pat_rec_max_outliers") );
  
  _inRecMode = true;
  
}

bool MINDfitman::in_rec_mode()
{
  return _inRecMode;
}

void MINDfitman::config_common_props()
{

  int vfit = _store.fetch_istore("vfit");
  int vnav = _store.fetch_istore("vnav");
  int vmod = _store.fetch_istore("vmod");
  int vmat = _store.fetch_istore("vmat");

  std::string info[4] = {"MUTE","NORMAL","VERBOSE","VVERBOSE"};

  l0 = Messenger::str(info[vfit]);
  l1 = Messenger::str(info[vnav]);
  l2 = Messenger::str(info[vmod]);
  l3 = Messenger::str(info[vmat]);

  _model = _store.fetch_sstore("model");
  _fitter = _store.fetch_sstore("kfitter");

}

void MINDfitman::config_geometry(Setup& det)
{

  _man.geometry_svc().set_zero_length( 0.1 * mm );
  _man.geometry_svc().set_infinite_length( 1e12 * mm );

  _man.geometry_svc().add_setup( "main", det );
  _man.geometry_svc().select_setup( "main" );

}


void MINDfitman::config_fitter()
{

  _man.fitting_svc()
    .select_fitter( _fitter );

  _man.fitting_svc()
    .set_fitting_representation( RP::slopes_curv_z );

  _man.fitting_svc()
    .fitter(_model).set_verbosity(l0);

}

void MINDfitman::config_navigator()
{

  _man.navigation_svc().navigator(_model)
    .set_unique_surface( true );
  
  _man.navigation_svc().set_verbosity(l1);
  _man.navigation_svc().navigator(_model)
    .set_verbosity(l1);

  _man.navigation_svc().inspector("ms")
    .set_verbosity(l1);

  _man.navigation_svc().navigator(_model)
    .master_inspector().set_verbosity(l1);

  _man.navigation_svc().inspector("BFieldMap")
    .set_verbosity(l1);

  _man.navigation_svc().inspector("eloss")
    .set_verbosity(l1);
  
  _man.navigation_svc().navigator(_model).set_max_number_steps(100);

  _man.navigation_svc().navigator(_model).set_allow_extrap_outside_volume(true);

  //_man.navigation_svc().navigator(_model).set_max_prop_distance_outside_finalvolume(150);


}

void MINDfitman::config_model()
{
  _man.model_svc().select_model(_model);

  _man.model_svc().model(_model)
    .equation().set_verbosity(l2);

  _man.model_svc().model(_model)
    .propagator().set_verbosity(l2);

  _man.model_svc().model(_model)
    .tool("noiser/ms").set_verbosity(l2);

  _man.model_svc().model(_model)
    .intersector("plane").set_verbosity(l1);

  _man.model_svc().model(_model)
    .intersector("numerical").set_verbosity(l1);
  
  _man.model_svc().model(_model)
    .intersector("helix_num").set_verbosity(l1);
  
  _man.model_svc().model(_model)
     .tool("correction/eloss").set_verbosity(l2);
}

void MINDfitman::config_matching()
{

  _man.matching_svc()
    .set_matching_representation( RP::slopes_curv_z );

  //Cellular automaton specifics
  _man.matching_svc()
    .retrieve_trajectory_finder<CellularAutomaton>("cat").set_verbosity(l3);


  _man.matching_svc().set_verbosity(l3);


  _man.matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_max_distance( _store.fetch_dstore("max_sep") * cm );

  _man.matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_plane_tolerance( sqrt(
			 pow(_store.fetch_dstore("pos_resX") * cm,2) +
			 pow(_store.fetch_dstore("pos_resY") * cm,2)));

  _man.matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_max_trajs( _store.fetch_istore("max_traj") );

}

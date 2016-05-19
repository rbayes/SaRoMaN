#include <MINDplotter.h>


class sortTrajByHits{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->size() > t2->size()) return true;
    return false;
  }

};
class sortTrajByLength{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->length() > t2->length()) return true;
    return false;
  }

};

//*************************************************************************************
MINDplotter::MINDplotter() {
  //*************************************************************************************

}

//*************************************************************************************
MINDplotter::~MINDplotter() {
  //*************************************************************************************
}

//*************************************************************************************
void MINDplotter::Initialize(string outFileName, bool patRec, bool clust,
			     bhep::prlevel vlevel) {
  //*************************************************************************************

  //bool ok = true;
  _patR = patRec;
  _clu = clust;

  _level = vlevel;

  _m = bhep::messenger(_level);

  _m.message("++Creating output root file++",bhep::VERBOSE);

  outFile = new TFile(outFileName.c_str(), "recreate");

  statTree = new TTree("tree", "Tree with pattern rec and fit data for MIND");

  statTree->SetAutoSave(500000);
  define_tree_branches();

  _m.message("++plotter initialized",bhep::VERBOSE);

  

  //return ok;
}

//*************************************************************************************
void MINDplotter::Execute(fitter& Fit, const bhep::event& evt) {
  //*************************************************************************************
  _m.message("++MINDplotter execute", bhep::VERBOSE);  
  bool ok1, ok2;


  ///event informations  
  _evNo = evt.event_number();

  ///vertex of the event
  bhep::Point3D evVert = evt.vertex();
  _vert[0] = evVert.x(); _vert[1] = evVert.y(); _vert[2] = evVert.z();
 

  if (_clu){
    get_tru_int_type( evt );
    if ( evt.find_iproperty("Charm") )
      _charm[0] = evt.fetch_iproperty("Charm");
    else _charm[0] = 0;
    if ( evt.find_iproperty("CharmHad") )
      _charm[1] = evt.fetch_iproperty("CharmHad");
    else _charm[1] = 0;
    if ( evt.find_dproperty("Q2") )
      _Q2 = evt.fetch_dproperty("Q2");
    else _Q2 = 0.;
    if ( evt.find_dproperty("EngTrans") )
      _engTrans = evt.fetch_dproperty("EngTrans");
    if ( evt.find_iproperty("nuType") )
      _pdg[0] = evt.fetch_iproperty("nuType");
    else _pdg[0] = 0;
    if ( evt.find_iproperty("intpart") )
      _pdg[1] = evt.fetch_iproperty("intpart");
    else _pdg[1] = 0;
    if ( evt.find_iproperty("nucType") )
      _pdg[2] = evt.fetch_iproperty("nucType");
    else _pdg[2] = 0;
  }



  
  ///Event Info
  _failEvent = 0;
  _failEvent = Fit.get_fail_event();

  // Initializing the  variable
  _trajs_no = 0; 
  _hitsInEvent = 0; 
  /// For hadrons
  _nhad[0] = 0; _nhad[1] = 0;
  // _tQ = 0;
  //   _tX[0]  = 0; _tX[1]  = 0;
  //   _tTh[0] = 0; _tTh[1] = 0;

  _hadE[0][0] = 0.; _hadE[1][0] = 0.; _hadE[2][0] = 0.; _hadE[3][0] = 0.;
  _hadE[0][1] = 0.; _hadE[1][1] = 0.; _hadE[2][1] = 0.; _hadE[3][1] = 0.;
  
  _hadP[0] = 0.;  _hadP[1] = 0.;  _hadP[2] = 0.;
  
  _hadEInNucleus = 0;
  
  _chadP[0] = 0.;  _chadP[1] = 0.;  _chadP[2] = 0.;  _chadP[3] = 0.;
  
  _nhadP[0] = 0.;  _nhadP[1] = 0.;  _nhadP[2] = 0.;  _nhadP[3] = 0.;

  
  _muontraj[0] = Fit.getMuonIndex()[0];  _muontraj[1] = Fit.getMuonIndex()[1];

  for (int i = 0;i<10;i++){
    _engTraj[i] = 0;
    _Fit[i]     = 0;
    _reFit[i]   = 0;
    _fail[i]    = 0;
    _vertZ[i]   = 0;
    _leng[i]    = 0;
    _ndf[i]     = 0;
    _nonMuonEdep[i] = 0;
    _hadHits[i] = 0; 
    _showerDir[i][0] = 0;   _showerDir[i][1] = 0;  _showerDir[i][2] = 0;


    _tQ[i] = 0; _tPDG[i] = 0;
    _tX[0][i]  = 0; _tX[1][i]  = 0;
    _tTh[0][i] = 0; _tTh[1][i] = 0;
    _tUnitDir[i][0] = 0; _tUnitDir[i][1] = 0; _tUnitDir[i][2] = 0;
    _recUnitDir[i][0] = 0; _recUnitDir[i][1] = 0; _recUnitDir[i][2] = 0;
   
    ///_hadE[1] = 0;
    
    _Xtent[i]=0;///
    _plns[0][i] = 0; _plns[1][i] = 0;
    _Q[0][i] = 0;  _Q[1][i] = 0;
    _rangqP[0][i] = 0;  _rangqP[1][i] = 0;
    _hitType[0][i] = 0;  _hitType[1][i] = 0;  _hitType[2][i] = 0;
    _qP[0][i] = 0; _qP[1][i] = 0;

    _positionStart[i][0]=0;
    _positionStart[i][1]=0;
    _positionStart[i][2]=0;
   
    
  }
  

  ///clear the vectors
  _XPos.clear();
  _YPos.clear();
  _ZPos.clear();
  _Edep.clear();
  _HTime.clear();
  _NodeFitted.clear();
  // _positionStart.clear();

  _XMeas.clear();
  _YMeas.clear();
  _ZMeas.clear();
  _EMeas.clear();
  _MuProp.clear();

  _XHadPos.clear();
  _YHadPos.clear();
  _ZHadPos.clear();
  

  for (int i = 0;i<2;i++){
    for (int j=0; j<2; j++){
      for (int k=0; k<10; k++){
	_X[i][j][k]=0;
	_Th[i][j][k]=0;
      }
    }
  }
  
  
  
    
  /// extract parameters for event only
  //get_event_info( evt,Fit);
  
  
 
  
  ///vector of trajectories
  std::vector<Trajectory*> &trajs = Fit.get_trajs();
  

  //total no of trajectories
  _trajs_no = trajs.size() < 10 ? trajs.size() : 10;

  ///length of trajectory
  double maxlength;
  int longesttraj;
  bool allfit=true;

  cout<<"In MINDplotter.cpp _clu="<<_clu<<endl;

  ///true particle informations
  if ( _clu )
    ok1 = extract_true_particle2(evt);
  // else
  //  ok1 = extract_true_particle1(evt, i);
  
  if( (int)_trajs_no == 0 ){
    int i = 0;
    // fill hadron information in the case of
    if(Fit.get_had_hits().size()==1){
      _hadHits[i] = Fit.get_had_hits()[i];
      // reconstructed unit direction for hadron shower
      _showerDir[i][0] = Fit.get_had_RecDir()[i][0];
      _showerDir[i][1] = Fit.get_had_RecDir()[i][1];
      _showerDir[i][2] = Fit.get_had_RecDir()[i][2];
      
      // edep for nonmuon hits
      _nonMuonEdep[i] = Fit.get_nonMuon_edep()[i];

      _vertZ[i] = Fit.showerVertZ(i);
      _Xtent[i] = Fit.showerXtent(i);
      _plns[0][i] = Fit.showerNplanes(i);
      _plns[1][i] = 0;
      _visEng = Fit.get_nonMuon_edep()[i];
    }
    else {
      _hadHits[0] = 0;
      _showerDir[0][0] = 0.0;
      _showerDir[0][1] = 0.0;
      _showerDir[0][2] = 1.0;
      _nonMuonEdep[0] = 0.0;
      _vertZ[0] = 0.0;
      _Xtent[0] = 0.0;
      _plns[0][0] = 0;
      _plns[1][0] = 0;
      _visEng = 0.0;
    }
  }

  //loop over trajectories
  for(int i=0; i<(int)_trajs_no; i++ ){
    
    _hadHits[i] = Fit.get_had_hits()[i];

    // reconstructed unit direction for hadron shower
    _showerDir[i][0] = Fit.get_had_RecDir()[i][0];
    _showerDir[i][1] = Fit.get_had_RecDir()[i][1];
    _showerDir[i][2] = Fit.get_had_RecDir()[i][2];
    
    // edep for nonmuon hits
    _nonMuonEdep[i] = Fit.get_nonMuon_edep()[i];

    //cout<<"inside MINDplotter:: Trajectory "<<i<<" and meas"<<*trajs[i]<<endl;

    ///fit and refit info  
    bool success  = trajs[i]->quality("fitted") && 
      (int)trajs[i]->quality("fitcheck") > 0;

    _Fit[i]       = (int)success;
    _fail[i]      = (int)(trajs[i]->quality("failType"));
    _reFit[i]     = (int)trajs[i]->quality("reseed");
    _vertZ[i]     = trajs[i]->quality("vertZ");
    _Xtent[i]     = (double)(trajs[i]->quality("xtent"));
    _ndf[i]       = trajs[i]->ndof();
    _m.message("_Fit =",_Fit[i]," _reFit =",_reFit[i],"  fail=",_fail[i]," vertexZ=  ",_vertZ[i],bhep::VERBOSE);    

    ///no of planes and freeplanes for each traj  
    
    if ( _failEvent != 7 ){
      _visEng = Fit.get_classifier().get_vis_eng();    
      
      _plns[0][i] = (int)(trajs[i]->quality("nplanes")) ;
      _plns[1][i] = (int)(trajs[i]->quality("freeplanes")) ;
      
      //cout<<"nplanes= "<<_plns[0][i]<<" freeplanes = "<<_plns[1][i]<<endl;
      
    } else { _plns[0][i] = 0; _plns[1][i] = -1; _visEng = 0; }
    
    
    
    ///interaction type
    if ( _failEvent != 7 )
      _intType[i] = (int)(trajs[i]->quality("intType")) ;
    else _intType[i] = 7;
    

    _rangqP[0][i] = 1./trajs[i]->quality("initialqP"); //Error probing.  

    ///for fitting success
    if (success) {
      State ste;       
      
      ///extract reconstructed params not extrapolating to vertex
      fill_kinematics(*trajs[i], ste, i);///
      
      ///extrapolate upto vertex

      ok2 = extrap_to_vertex(*trajs[i], evt.vertex(), Fit, ste, i);

      if(trajs[i]->quality("lowPt") == 1)  // Get around not having a vertex.
	{
	  ok2 = true;
	}
      // 
      
      ///if extrapolate to vertex is successful then fill the parameters
      if (ok2) {
	position_pulls(i);
	direction_pulls(i);
	momentum_pulls(i);
      }
      
      ///length calculation :: No need to change the "sign" because after refiting
      /// the nodes are reversed back inside event classifier
      _leng[i] = trajs[i]->length();

      ///longest track
      if(abs(_leng[i])>maxlength && _leng[i] > 50){
	maxlength = _leng[i];
	longesttraj = i;
      } else if(_leng[i] < 50)
	allfit = false;
      
      ///range of q/p calculation ?

      _rangqP[0][i] = 1./trajs[i]->quality("initialqP"); //Original
      // _rangqP[1][i] = Fit.get_classifier().RangeMomentum(trajs[i]->length());
      // _rangqP[1][i] = _rangqP[0][i] - _qP[0][i];
      // _rangqP[2][i] = _qP[1][i] - _rangqP[0][i];
      
      ///chi2 value of traj
      max_local_chi2( *trajs[i], i );
      
    
    }
  
    //If fit not successful set rec values to zero.
    if (!success){
      allfit=false;
      _X[0][0][i] = 0; _X[0][1][i] = 0; 
      _X[1][0][i] = 0; _X[1][1][i] = 0; 
      //_X[3][0][i] = 0; _X[3][1][i] = 0; 
      //_X[4][0][i] = 0; _X[4][1][i] = 0; 
      //_X[5][0][i] = 0; _X[5][1][i] = 0;

      _Th[0][0][i] = 0; _Th[0][1][i] = 0;
      _Th[1][0][i] = 0; _Th[1][1][i] = 0;
      //_Th[3][0][i] = 0; _Th[3][1][i] = 0;
      //_Th[4][0][i] = 0; _Th[4][1][i] = 0;
      //_Th[5][0][i] = 0; _Th[5][1][i] = 0;

      _qP[0][i] = 0; _qP[1][i] = 0; // _qP[3][i] = 0; _qP[4][i] = 0; _qP[5][i] = 0;

      _Chi[0][i] = 0; _Chi[1][i] = 0; _Chi[2][i] = 0;

      _Q[0][i] = 0; 
      _Q[1][i] = 0;

    }

    
  }//end of traj loop

 /// reconstruct the hadron direction
  int muonsel = _muontraj[0];
  if(allfit) muonsel = longesttraj;
  // hadron_direction(Fit);

  /// 
  if (_patR){
    if ( _clu ) {
      hitBreakUp(Fit);
      patternStats2( Fit);
      //  else
      // patternStats1( Fit, *trajs[i],i );
    }
  }
  
  

  //Fill tree event with the values.
  //int fillcatch = statTree->Fill();
  statTree->Fill();


  // clear vector of trajectories
  trajs.clear();
 

  // if(evt % 1000) statTree->AutoSave("FlushBaskets");
  //return ok1;
}


/*****************************************************************/
void MINDplotter::get_tru_int_type(const bhep::event& evt){
  
  string intName = evt.fetch_sproperty("IntType");
  
  if ( intName=="CCQE" || intName=="<QES - Weak[CC]>" )
    _truInt = 1;
  else if ( intName=="NCQE" || intName=="<QES - Weak[NC]>" )
    _truInt = 2;
  else if ( intName=="CCDIS" || intName=="<DIS - Weak[CC]>" )
    _truInt = 3;
  else if ( intName=="NCDIS" || intName=="<DIS - Weak[NC]>" )
    _truInt = 4;
  else if ( intName=="1piRes" )
    _truInt = 5;
  else if ( intName=="miscRes" || intName=="<RES - Weak[CC]>" || intName=="<RES - Weak[NC]>" )
    _truInt = 6;
  else if ( intName=="eEl" || intName=="<NuEEL - Weak[NC]>" || intName=="<NuEEL - Weak[CC]>" )
    _truInt = 7;
  else if ( intName=="muINVe" || intName=="<IMD - Weak[CC]>" )
    _truInt = 7;
  else
    _truInt = 8;
  
}

//*************************************************************************************
void MINDplotter::Finalize() {
  //*************************************************************************************

  //bool ok = true;

  _m.message("++Finalizing Output++",bhep::VERBOSE);
  
  outFile->Write();
  outFile->Close();

  delete outFile;
  // if(_intType==5){
  // _CAFile->Write();
  // _CAFile->Close();//}

  //return ok;
}

//*************************************************************************************
void MINDplotter::define_tree_branches() {
  //*************************************************************************************

  ///for each event
  statTree->Branch("MC_Evt", &_evNo, "MC_EventNo/I");
  statTree->Branch("MC_evVertex", &_vert, "MC_vert[3]/D");
  statTree->Branch("recon_EvtFail", &_failEvent, "recon_EvtFail/I");
  if (_clu) statTree->Branch("MC_TrueInteraction",&_truInt,"MC_truInt/I");
  statTree->Branch("MC_PrimaryEng", &_nuEng, "MC_KEng/D");
  if (_clu) {
    //statTree->Branch("Charm", &_charm, "isCharm/I:pdg/I");
    //statTree->Branch("InteractionQ2", &_Q2, "Q2/D");
    //statTree->Branch("EventPdg",&_pdg, "initnu/I:partPDG/I:nucType/I");
  }
  statTree->Branch("raw_visibleEng", &_visEng, "raw_visEng/D");
  statTree->Branch("raw_HitsInEvent", &_hitsInEvent, "raw_TotalHits/I");
 
  ///for each traj 
  statTree->Branch("classif_NumTrajectory", &_trajs_no, "trajNo/I");
  statTree->Branch("traj_hadronHits", &_hadHits, "traj_hadHits[trajNo]/I");
  statTree->Branch("traj_NonMuonEdep", &_nonMuonEdep, "traj_HadEdep[trajNo]/D");
  statTree->Branch("classif_LongMuTraj", &_muontraj[0], "longMuTraj/I");
  statTree->Branch("classif_NoMuTraj", &_muontraj[1], "noMuTraj/I");

  statTree->Branch("traj_TrajVertex", &_vertZ, "traj_trajVert[trajNo]/D");
  statTree->Branch("traj_Fitted", &_Fit, "traj_success[trajNo]/I");
  statTree->Branch("traj_backFit",&_reFit,"traj_backFit[trajNo]/I");
  statTree->Branch("traj_Fail", &_fail, "traj_FailType[trajNo]/I");
  statTree->Branch("traj_interaction",&_intType,"traj_Inter[trajNo]/I");

  ///
  statTree->Branch("traj_visEng",&_engTraj, "traj_eng[trajNo]/D");
  // statTree->Branch("visEngSel", &_engSelTraj, "engSelTraj[trajNo]/D");
  statTree->Branch("traj_EngDep",&_engvar[0],"traj_meanDep[trajNo]/D");
  statTree->Branch("traj_EngVar",&_engvar[1],"traj_engVar[trajNo]/D");
  statTree->Branch("traj_LowDep",&_engvar[2],"traj_lowDep[trajNo]/D");
  statTree->Branch("traj_HighDep",&_engvar[3],"traj_highDep[trajNo]/D");

  
  statTree->Branch("traj_Planes", &_plns[0], "traj_nplanes[trajNo]/I");
  statTree->Branch("traj_FPlanes",&_plns[1], "traj_freeplanes[trajNo]/I");
 

  ///position  
  statTree->Branch("MCtr_TrueXPosition", &_tX[0][0], "MCtr_vertXPos[trajNo]/D");
  statTree->Branch("MCtr_TrueYPosition", &_tX[1][0], "MCtr_vertYPos[trajNo]/D");
  statTree->Branch("traj_vertXPosition", &_X[0][0],"traj_vertXPos[trajNo]/D");
  statTree->Branch("traj_vertYPosition", &_X[0][1], "traj_vertYPos[trajNo]/D");
  statTree->Branch("traj_ErrXPosition", &_X[1][0], "traj_ErrXPos[trajNo]/D");
  statTree->Branch("traj_ErrYPosition", &_X[1][1], "traj_ErrYPos[trajNo]/D");
  //statTree->Branch("RecXPosition_WVE", &_X[3][0], "recXPos_WVE[trajNo]/D");
  //statTree->Branch("RecYPosition_WVE", &_X[3][1], "recYPos_WVE[trajNo]/D");
  //statTree->Branch("ErrXPosition_WVE", &_X[2][0], "ErrXPos_WVE[trajNo]/D");
  //statTree->Branch("ErrYPosition_WVE", &_X[2][1], "ErrYPos_WVE[trajNo]/D");
  //statTree->Branch("XPull", &_X[5][0], "pull_x[trajNo]/D");
  //statTree->Branch("YPull", &_X[5][1], "pull_y[trajNo]/D");

  ///Directions
  statTree->Branch("MCtr_XDirection", &_tTh[0], "MCtr_XTh[trajNo]/D");
  statTree->Branch("MCtr_YDirection", &_tTh[1], "MCtr_YTh[trajNo]/D");
  statTree->Branch("MCtr_Direction", _tUnitDir, "MCtr_Dir[trajNo][3]/D");
  statTree->Branch("traj_Direction", _recUnitDir, "traj_Dir[trajNo][3]/D");
  statTree->Branch("traj_XDirection", &_Th[0][0], "traj_XTh[trajNo]/D");
  statTree->Branch("traj_YDirection", &_Th[0][1], "traj_YTh[trajNo]/D");
  statTree->Branch("traj_ErrXDirection", &_Th[1][0], "traj_ErrXTh[trajNo]/D");
  statTree->Branch("traj_ErrYDirection", &_Th[1][1], "traj_ErrYTh[trajNo]/D");
  // statTree->Branch("RecXDirection_WVE", &_Th[3][0], "recXTh_WVE[trajNo]/D");
  // statTree->Branch("RecYDirection_WVE", &_Th[3][1], "recYTh_WVE[trajNo]/D");
  // statTree->Branch("ErrXDirection_WVE", &_Th[4][0], "ErrXTh_WVE[trajNo]/D");
  // statTree->Branch("ErrYDirection_WVE", &_Th[4][1], "ErrYTh_WVE[trajNo]/D");
  // statTree->Branch("XPullTh", &_Th[5][0], "pull_Thx[trajNo]/D");
  // statTree->Branch("YPullTh", &_Th[5][1], "pull_Thy[trajNo]/D");


  ///Momentum
  statTree->Branch("MCtr_Mom", &_tqP, "MCtr_qP[trajNo]/D");
  statTree->Branch("traj_Mom", &_qP[0], "traj_qP[trajNo]/D");
  statTree->Branch("traj_ErrMom", &_qP[1], "traj_ErrqP[trajNo]/D");
  statTree->Branch("traj_MomNE", &_qP[2], "traj_qPNE[trajNo]/D");
  // statTree->Branch("RecMom_WVE", &_qP[3], "recqP_WVE[trajNo]/D");
  // statTree->Branch("ErrMom_WVE", &_qP[4], "ErrqP_WVE[trajNo]/D");
  // statTree->Branch("MomPull", &_qP[5], "pull_mom[trajNo]/D");



  ///Charge 
  statTree->Branch("MCtr_Charge", &_tQ, "MCtr_Q[trajNo]/I");
  statTree->Branch("MCtr_PDG", &_tPDG, "MCtr_truPDG[trajNo]/I");
  statTree->Branch("traj_Charge", &_Q[0], "traj_Q[trajNo]/I");
  statTree->Branch("MCtr_ID", &_Q[1], "MCtr_ID[trajNo]/I");



  statTree->Branch("traj_length", &_leng,"traj_len[trajNo]/D");
  ///range of qP calculations only for success ?
  statTree->Branch("traj_initrangqP", &_rangqP[0],"traj_intrangqP[trajNo]/D");
  // statTree->Branch("traj_rangqP", &_rangqP[1],"traj_rangqP[trajNo]/D");
  //statTree->Branch("rangErr", &_rangqP[1],"rangErr[trajNo]/D");
  //statTree->Branch("rangdiff", &_rangqP[2],"recrangediff[trajNo]/D");

  statTree->Branch("traj_FitChiInfo", &_Chi[0], "traj_Chi[trajNo]/D");
  statTree->Branch("traj_MaxFitChiInfo", &_Chi[1], "traj_MaxLoc[trajNo]/D");
  statTree->Branch("traj_TotalFitChiInfo", &_Chi[2], "traj_TotalChi2[trajNo]/D");
  statTree->Branch("traj_NDoF", &_ndf, "traj_ndof[trajNo]/I");
  statTree->Branch("traj_extent", &_Xtent, "traj_xtent[trajNo]/D");

  
  statTree->Branch("traj_NoHits", &_nhits, "traj_nhits[trajNo]/I");
  statTree->Branch("classif_Nallhits", &_nallhits, "classif_traj_nallhits/I"); 
  statTree->Branch("trajNode_XPos", &_XPos,32000,0);
  statTree->Branch("trajNode_YPos", &_YPos, 32000,0);
  statTree->Branch("trajNode_ZPos", &_ZPos,32000,0);
  statTree->Branch("trajNode_EngDeposit", &_Edep,32000,0);
  statTree->Branch("trajNode_Time", &_HTime,32000,0);
  statTree->Branch("trajNode_Fitted", &_NodeFitted,32000,0);
  statTree->Branch("traj_positionStart",_positionStart,"traj_positionStart[trajNo][3]/D" );
  
  statTree->Branch("raw_Xmeas", &_XMeas,32000,0);  
  statTree->Branch("raw_Ymeas", &_YMeas, 32000,0);
  statTree->Branch("raw_Zmeas", &_ZMeas,32000,0);
  statTree->Branch("raw_EngMeas", &_EMeas,32000,0);
  statTree->Branch("raw_MotherProp", &_MuProp, 32000, 0);

  statTree->Branch("raw_MotherMeas", &_MotherMeas,32000,0);
  statTree->Branch("raw_MotherProp", &_MotherProp,32000,0);

  //statTree->Branch("theta0", &_Theta0,32000,0);
  statTree->Branch("raw_theta1", &_Theta1,32000,0);
  statTree->Branch("raw_theta2", &_Theta2, 32000,0);
  statTree->Branch("raw_theta3", &_Theta3, 32000,0);

  statTree->Branch("raw_delta1", &_Delta1,32000,0);
  statTree->Branch("raw_delta2", &_Delta2, 32000,0);

  ///for hadrons
  //statTree->Branch("NChadnNhad", &_nhad[0], "nChad/I");
  //statTree->Branch("Nhad", &_nhad[1], "nNhad/I");
  //statTree->Branch("hadronP", &_had, "hadP[3]/D");
  //statTree->Branch("EnergyTransfer", &_engTrans, "nu/D");
  //statTree->Branch("hadTruEng", &_hadE[0][0], "truE/D");
  //statTree->Branch("hadTrueE", &_hadEInNucleus, "hadTrueE/D");

  //statTree->Branch("hadRdirP", &_hrecP, "hrdirP[2][3]/D");
  // statTree->Branch("hadRecEng", &_hadE[1], "hrecE[2]/D");
  //statTree->Branch("hadQESEng", &_hadE[2], "hQESE[2]/D");
  // statTree->Branch("hadtransPrfl", &_hadE[3], "htprfl[2]");
  //statTree->Branch("RecShowerDir", _showerDir, "showerDir[10][3]/D");

  //statTree->Branch("ChadMomX", &_chadP[0], "chadPx/D");
  //statTree->Branch("ChadMomY", &_chadP[1], "chadPy/D");
  //statTree->Branch("ChadMomZ", &_chadP[2], "chadPz/D");
  //statTree->Branch("ChadMomE", &_chadP[3], "chadE/D");

  //statTree->Branch("NhadMomX", &_nhadP[0], "nhadPx/D");
  //statTree->Branch("NhadMomY", &_nhadP[1], "nhadPy/D");
  //statTree->Branch("NhadMomZ", &_nhadP[2], "nhadPz/D");
  //statTree->Branch("NhadMomE", &_nhadP[3], "nhadE/D");
  
  //statTree->Branch("hadDir", &_haddot, "dotProd[2]/D");
  //statTree->Branch("truHitIndex", &_truHitIndex[0], "truHitInd[trajNo]/I");
  //statTree->Branch("candHitIndex", &_truHitIndex[1], "candHitInd[trajNo]/I");
  

  ///HitBreakUp
  ///statTree->Branch("TruMu", &_hitType[0], "nTruMu/I");
  statTree->Branch("MC_TruMu", &_hitTrMu, "MC_nTruMu/I");
  statTree->Branch("traj_InMu", &_hitType[0], "traj_nInMu[trajNo]/I");
  statTree->Branch("MCtr_MuInMu", &_hitType[1], "MCtr_nMuInMu[trajNo]/I");
  statTree->Branch("traj_fittedNode", &_hitType[2], "traj_nFitN[trajNo]/I");
  statTree->Branch("classif_hadN", &_hitHad, "classif_nhad/I");
  
  
  //statTree->Branch("MuHits", &_mus, "truMu[trajNo][nallhits]/I");
  //statTree->Branch("CandHits", &_cand, "inMu[trajNo][nallhits]/B");
  //statTree->Branch("FittedNodes",&_node,"fitNode[trajNo][nallhits]/B");
  //statTree->Branch("HadronNodes",&_had,"hadN[trajNo][nallhits]/B");
  statTree->Branch("traj_PatRecChi", &_pChi, "maxChiMu[trajNo]/D:MinChiHad[trajNo]/D:MaxConsecHol[trajNo]/D");


  statTree->Branch("classif_XHadPos", &_XHadPos,32000,0);
  statTree->Branch("classif_YHadPos", &_YHadPos, 32000,0);
  statTree->Branch("classif_ZHadPos", &_ZHadPos,32000,0);
  
  
}

//*************************************************************************************
void MINDplotter::fill_kinematics(const Trajectory& traj, State& ste, const int trajno){ 
  //*************************************************************************************

  _m.message("++Fill vertex kinematics++",bhep::VERBOSE);

  /// 1st fitted state of the traj
  ste = traj.node(traj.first_fitted_node()).state();

  // information of the 1st fitted node
  EVector v = ste.hv().vector();
  EMatrix C = ste.hv().matrix();

  //Reconstructed x,y position.
  //_X[3][0][trajno] = v[0]; _X[3][1][trajno] = v[1];
 

  //Corresponding Error.
  //if (C[0][0]>0)
  //  _X[4][0][trajno] = sqrt(C[0][0]);
  //if (C[1][1]>0)
  //  _X[4][1][trajno] = sqrt(C[1][1]);
  
 
  //Reconstructed q/P without extrapolating to vertex
  if (v[5] !=0) _Q[0][trajno] = (int)( v[5]/fabs(v[5]) );
  _qP[2][trajno] = 1./v[5];
  
  //cout<<"recqP_WVE="<< _qP[2][trajno]<<endl;
  //Corresponding Error.
  //if (C[5][5]>0){
  //  _qP[4][trajno] = sqrt(C[5][5]);
  
  //momentum pull
  //_qP[5][trajno] = (_qP[3][trajno] - _qP[0][trajno])/_qP[4][trajno];
  // }
  
  //Correctly ID'd charge?.
  if (_tQ[trajno] == _Q[0][trajno]) _Q[1][trajno] = 1;
  else _Q[1][trajno] = 0;
  /*
  //cout<<" fill_kinematics: truQ="<<_Q[0][trajno]<<"  recQ="<<_Q[1][trajno]<<"  ID="<<_Q[2][trajno]<<endl;
  //Reconstructed direction.
  _Th[3][0][trajno] = v[3]; _Th[3][1][trajno] = v[4];
  
  //Corresponding error.
  if (C[3][3]>0){
    _Th[4][0][trajno] = sqrt(C[3][3]);
    //direction pull
    _Th[5][0][trajno] = (_Th[3][0][trajno] - _Th[0][0][trajno])/_Th[4][0][trajno] ;}
  
  
  if (C[4][4]>0){
    _Th[4][1][trajno] = sqrt(C[4][4]);
    //direction pull
    _Th[5][1][trajno] = (_Th[3][1][trajno] - _Th[0][1][trajno])/_Th[4][1][trajno] ;
  }
  */
  
}

//*************************************************************************************
bool MINDplotter::extrap_to_vertex(const Trajectory& traj, 
				   const bhep::Point3D& vertexLoc,
				   fitter& fitObj, State& ste, const int trajno) {
  //*************************************************************************************

  _m.message("++Extrapolation function, Finding best fit to vertex++",bhep::VERBOSE);
  

  ///
  _fittedVert = EVector(3,0);
  _vertMat = EMatrix(3,3,0);
  // ste.clear();

  ///set the 1st fitted node state
  ste = traj.node(traj.first_fitted_node()).state();

  /// start position for each trajectory
  for(int i=0; i<3; i++) {
    _positionStart[trajno][i]= ste.hv().vector()[i];
  }
   _m.message("start position =",_positionStart[trajno][0],",",_positionStart[trajno][1],",",_positionStart[trajno][2],bhep::VERBOSE);

  /// create the surface of the vertex
  EVector pos(3,0); pos[2] = vertexLoc.z();
  EVector axis(3,0); axis[2] = 1;

  double R1 = 1000000000;
  double R2 = 0;
  double l;

  Ring surf(pos, axis, R1, R2);

  //cout<<"1st fitted node Z="<<ste.vector()[2]<<"  ;tru vextexZ="<<pos[2]<<endl;

  bool ok;

  if(traj.quality("lowPt") != 1) // Get around not having a vertex.
    {
      /// Add the surfaceof vertex and prpagate to that surface
      fitObj.man().geometry_svc().setup().add_surface("mother","vertex",&surf);
      ok = fitObj.man().navigation_svc().propagate(surf,ste,l);
      fitObj.man().geometry_svc().setup().remove_surface("vertex");
      
      //Convert to slopes representation.
      // fitObj.man().model_svc().conversion_svc().representation().convert(ste, RP::slopes_curv_z);
    }
  
  //Grab fitted vertex information.
  _fittedVert = ste.hv().vector();
  _vertMat = ste.hv().matrix();
  _m.message("extrapolated vertex =",_fittedVert[2],"\n",bhep::DETAILED);

  return ok;
}


//**************************************************************************************
void MINDplotter::position_pulls(const int trajno) {
  //**************************************************************************************

  //Function to calculate the pull for a measurement
  _m.message("++Calculating position pulls++",bhep::VERBOSE);
  
  //Reconstructed x,y position.
  _X[0][0][trajno] = _fittedVert[0]; _X[0][1][trajno] = _fittedVert[1];
  
  //cout<<"recPos ="<<_X[1][0][trajno]<<endl;
  
  
  //Corresponding Error.
  if (_vertMat[0][0]>0){
    _X[1][0][trajno] = sqrt(_vertMat[0][0]);
    
    
    //pull of X
    // _X[5][0][trajno] = (_X[1][0][trajno] - _X[0][0][trajno])/_X[2][0][trajno] ;
  }
  if (_vertMat[1][1]>0){
    _X[1][1][trajno] = sqrt(_vertMat[1][1]);
    
    
    //pull of Y
    // _X[5][1][trajno] = (_X[1][1][trajno] - _X[0][1][trajno])/_X[2][1][trajno] ;
  }
  
  // cout<<"rec position ="<< _X[0][0][trajno]<<"   "<< _X[0][1][trajno]<<endl;
  
}

//**************************************************************************************
void MINDplotter::momentum_pulls(const int trajno) {
  //**************************************************************************************

  ///Function to calculate momentum pulls.
  _m.message("++Calculating momentum pulls++",bhep::VERBOSE);

  //Reconstructed q/P.
  if (_fittedVert[5] !=0) _Q[0][trajno] = (int)( _fittedVert[5]/fabs(_fittedVert[5]) );
  _qP[0][trajno] = 1./_fittedVert[5];

  _m.message("rec Mom =",1./_fittedVert[5],bhep::VERBOSE);
  
  //Corresponding Error.
  if (_vertMat[5][5]>0)
    _qP[1][trajno] = sqrt(_vertMat[5][5])/_fittedVert[5]/_fittedVert[5];
  
  //Correctly ID'd charge?.
  if (_tQ[trajno] == _Q[0][trajno]) _Q[1][trajno] = 1;
  else _Q[1][trajno] = 0;
  
}

//**************************************************************************************
void MINDplotter::direction_pulls(const int trajno) {
  //**************************************************************************************

  //Reconstructed direction.
  _Th[0][0][trajno] = _fittedVert[3]; _Th[0][1][trajno] = _fittedVert[4];


  //Reconstructed unit direction
  double mag = sqrt(1 + pow(_fittedVert[3],2) + pow(_fittedVert[4],2));
  _recUnitDir[trajno][0] = _fittedVert[3]/mag;
  _recUnitDir[trajno][1] = _fittedVert[4]/mag;
  _recUnitDir[trajno][2] = 1./mag;
  

  //Corresponding error.
  if (_vertMat[3][3]>0)
    _Th[1][0][trajno] = sqrt(_vertMat[3][3]);
 
  if (_vertMat[4][4]>0)
    _Th[1][1][trajno] = sqrt(_vertMat[4][4]);

}

/*//*************************************************************************************
void MINDplotter::get_event_info(const bhep::event& evt,fitter& fit){
//*************************************************************************************
  _hitTrMu = 0;
 
  
  /// neutrino energy of the particle
  _nuEng = evt.fetch_dproperty("nuEnergy") * bhep::MeV;
  
  /// measurement vector
  std::vector<cluster*>& meas = fit.GetMeasVec();
  
  //cout<<" meas= "<<fit.GetMeas(0)->position()[2]<<"  "<<fit.GetNMeas()<<" vector size= "<<meas.size()<<endl;
  
  
  ///loop over measurements to calculate toltal no of truMu
  for(int i=0; i<fit.GetNMeas(); i++ ){
    // cout<<"event_info::  pos="<<fit.GetMeas(i)->position()[2]<<endl;
    
    if ( meas[i]->get_mu_prop() > 0.8 )  _hitTrMu++;
    
  }
  _m.message("Total TruMu number =",_hitTrMu, bhep::VERBOSE);
}

*/
  

//*************************************************************************************
bool MINDplotter::extract_true_particle1(const bhep::event& evt, const int trajNo) {

  //*************************************************************************************
  
  /* sets true particle momentum for calculation and returns a reference
     to the particle */

  //_nuEng = atof( evt.fetch_property("Enu").c_str() ) * bhep::GeV;
  //_nuEng = evt.fetch_dproperty("nuEnergy") * bhep::MeV;
  
  ///vector of true particles
  const vector<bhep::particle*> Pospart = evt.true_particles();
  
  /// loop over true particles
  int count = 0;
  for (int iParts=0;iParts < (int)Pospart.size();iParts++){
    if (Pospart[iParts]->name().compare("mu-")==0 ||
	Pospart[iParts]->name().compare("pi-")==0 ||
	Pospart[iParts]->name().compare("e-")==0){
      _tQ[trajNo] = -1;
      _truPart = Pospart[iParts];
      count++;
    } 
    else if (Pospart[iParts]->name().compare("mu+")==0 ||
	     Pospart[iParts]->name().compare("pi+")==0 ||
	     Pospart[iParts]->name().compare("e+")==0){
      _tQ[trajNo] = 1;
      _truPart = Pospart[iParts];
      count++;
    } 
    if (Pospart[iParts]->name().compare("Hadronic_vector")==0){
      _hadP[0] = Pospart[iParts]->px();
      _hadP[1] = Pospart[iParts]->py();
      _hadP[2] = Pospart[iParts]->pz();
      _hadE[0][0] = atof( Pospart[iParts]->fetch_property("HadE").c_str() ) * bhep::GeV;
      // _hadP[3] = Pospart[iParts]->fetch_dproperty("HadE") * bhep::GeV;
    }
  }
  //STUFF ABOVE FOR REDESIGN.!!!!!

  if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
    _tQ[trajNo]  = 0;
    _tqP[trajNo] = 0;
    _tTh[0][trajNo] = _tTh[1][trajNo]= 0;
    _tX[0][trajNo]  = _tX[1][trajNo] =0;///
    return false;
  }
  
  //Set true values of muon mom. etc.
  //True q/P.
  _tqP[trajNo] = _tQ[trajNo]*_truPart->p();
  //True direction.
  _tTh[0][trajNo]= _truPart->px()/_truPart->pz();
  _tTh[1][trajNo]= _truPart->py()/_truPart->pz();
  
  //True x,y position.
  _tX[0][trajNo] = evt.vertex().x(); _tX[1][trajNo] = evt.vertex().y();
  
  // Pospart.clear();
  return true;
}




//**********************************************************************************
bool MINDplotter::extract_true_particle2(const bhep::event& evt) {
  //*************************************************************************************
 

  /// neutrino energy of the event
  _nuEng = evt.fetch_dproperty("nuEnergy") * bhep::MeV;
 
  _m.message(" Neurino Energy=",evt.fetch_dproperty("nuEnergy"), bhep::VERBOSE);

  ///vector of true particles
  const vector<bhep::particle*> Pospart = evt.true_particles();
  
  /*//----------------------------
  std::vector<int> vpdg = evt.fetch_ivproperty("fspdg");
  for(int np=0; np<(int)vpdg.size(); np++){
    cout<<"pdg="<<vpdg[np]<<endl;
  }
  vpdg.clear();
  //----------------------------*/


  ///loop over the true particles to extract info
  int count = 0, trackNo=0, had_count=0;

  cout<<"Pospart.size()="<<(int)Pospart.size()<<endl;

  for (int iParts=(int)Pospart.size()-1;iParts >= 0;iParts--){
    count=0;

    cout<<"in hadron="<<Pospart[iParts]->p()<<"  pdg="<<Pospart[iParts]->pdg()<<endl;
    //cout<<"ID="<<Pospart[iParts]->fetch_sproperty("G4TrackID")<<" parentID="<<Pospart[iParts]->find_dproperty("G4ParentID")<<endl;
	
   
    if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
	 Pospart[iParts]->find_sproperty("PrimaryLepton") ){
      _m.message("For primary lepton, mom =",Pospart[iParts]->p(),"  pdg=",Pospart[iParts]->pdg(), bhep::VERBOSE);
     
      cout<<"For primary lepton, mom ="<<Pospart[iParts]->p()<<"  pdg="<<Pospart[iParts]->pdg()<<endl;

      _truPart = Pospart[iParts];
	
      count++;
      trackNo = count;
    } else if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
		!Pospart[iParts]->find_sproperty("PrimaryLepton") ) 
      {
	///to check whether same particle or not
	if(Pospart[iParts] == Pospart[iParts+1]) had_count++;
	
	///1. energy >1 GeV (rough assumption), 2. if same particle has min 5 hits & 3. length of track (atleast 5 layers each 52 mm thick)
	if(Pospart[iParts]->p()>= 1*bhep::GeV || had_count>=5 || Pospart[iParts]->fetch_dproperty("length") >= 26*bhep::cm)  {
	  count ++;
	  trackNo = trackNo + count;
	  _truPart = Pospart[iParts];
	  _m.message("For hadron=",Pospart[iParts]->p(),"  pdg=",Pospart[iParts]->pdg()," trackNo=",trackNo, bhep::VERBOSE);
	}
	add_to_hads( *Pospart[iParts] );///
	cout<<"in hadron="<<Pospart[iParts]->p()<<"  pdg="<<Pospart[iParts]->pdg()<<endl;
	cout<<"ID="<<Pospart[iParts]->fetch_sproperty("G4TrackID")<<" parentID="<<Pospart[iParts]->find_dproperty("G4ParentID")<<endl;
	
      }

    if(count!=0) {
      
      int trajNo = trackNo -1;
      //charge
      _tQ[trajNo] = _truPart->charge() ;
      
       //pdg
      _tPDG[trajNo] = _truPart->pdg() ;

      //Set true values of muon mom. etc.
      //True q/P.
      _tqP[trajNo] = _tQ[trajNo]*_truPart->p();
    
      // cout<<"tru mom="<<1./_tqP[trajNo]<<" charge="<<_tQ[trajNo]<<endl;
      
      //True direction.
      _tTh[0][trajNo]= _truPart->px()/_truPart->pz();
      _tTh[1][trajNo]= _truPart->py()/_truPart->pz();
      
      //True Unit direction.
      double norm = sqrt(1 + pow( _tTh[0][trajNo],2) + pow(_tTh[1][trajNo],2));
      _tUnitDir[trajNo][0]=  _tTh[0][trajNo]/norm;
      _tUnitDir[trajNo][1]=  _tTh[1][trajNo]/norm;
      _tUnitDir[trajNo][2]=  1./norm;
      
      //True x,y position.
      _tX[0][trajNo] = evt.vertex().x(); 
      _tX[1][trajNo] = evt.vertex().y();
      
      _m.message("_tQ[",trajNo,"]=",_tQ[trajNo],"   ;tru mom=",_truPart->p()," 1/_tqP[",trajNo,"]=",1./_tqP[trajNo], bhep::VERBOSE);
      
      /*if ( Pospart[iParts]->fetch_dproperty("length") >= 30*bhep::cm) {
      trackNo++;
      cout<<"mom="<<Pospart[iParts]->p()<<endl;
      }*/
    } 
  }
  
  // cout<<"track no="<<trackNo<<endl;
  
  std::vector<double> hadInf = evt.fetch_dvproperty("had4vec");
  
  for (int itn = 0;itn < 4;itn++)
    _hadP[itn] = hadInf[itn]*bhep::GeV;
  _hadE[0][0] = hadInf[3]*bhep::GeV;
  hadInf.clear();
  

  //true hadron energy which are inside nucleus
  _hadEInNucleus = evt.fetch_dproperty("hadEInNucleus");
  
  
  /*if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
    _tQ[count]     = 0;
    _tqP[count]    = 0;
    _tTh[0][count] = _tTh[1][count]= 0;
    _tX[0][count]  = _tX[1][count] = 0;///
    
     
    return false;
    }*/
  
  
  
  /*
//-----------------------------------------------------
///loop over the true particles to extract info
  int count = 0;
  for (int iParts=(int)Pospart.size()-1;iParts >= 0;iParts--){
   
    if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
	 Pospart[iParts]->find_sproperty("PrimaryLepton") ){

      //cout<<"check="<< Pospart[iParts]->name()<<endl;
      if ( Pospart[iParts]->name() == "mu+"){
	_truPart = Pospart[iParts];
	_tQ[count] = 1;
	count++;
      } else if ( Pospart[iParts]->name() == "mu-"){
	_truPart = Pospart[iParts];
	_tQ[count] = -1;
	count++;
      }
    } else {
      if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
	   !Pospart[iParts]->find_sproperty("PrimaryLepton") ) 
	
	_truPart = Pospart[iParts];
      add_to_hads( *Pospart[iParts] );///
      cout<<"in hadron"<<endl;
     
    }

    std::vector<double> hadInf = evt.fetch_dvproperty("had4vec");
  
  for (int itn = 0;itn < 4;itn++)
    _hadP[itn] = hadInf[itn]*bhep::GeV;
  _hadE[0][0] = hadInf[3]*bhep::GeV;
  hadInf.clear();
  
  
  if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
    _tQ[trajNo]     = 0;
    _tqP[trajNo]    = 0;
    _tTh[0][trajNo] = _tTh[1][trajNo]= 0;
    _tX[0][trajNo]  = _tX[1][trajNo] = 0;///
    
     
    return false;
  }
  //Set true values of muon mom. etc.

  //True q/P.
  _tqP[trajNo] = _tQ[trajNo]/_truPart->p();

  cout<<"tru mom="<<1./_tqP[trajNo]<<endl;
  
  //True direction.
  _tTh[0][trajNo]= _truPart->px()/_truPart->pz();
  _tTh[1][trajNo]= _truPart->py()/_truPart->pz();
  
  //True x,y position.
  _tX[0][trajNo] = evt.vertex().x(); 
  _tX[1][trajNo] = evt.vertex().y();*/
  //--------------------------------------------------------------
 
  
  // Pospart.clear();
  // cout<<"tru position ="<< _X[0][0][trajno]<<"   "<< _X[0][1][trajno]<<endl;
  // std::cout<<"end of extract_true_particle"<<std::endl;
  return true;
}



/**************************************************************/
void MINDplotter::add_to_hads(const bhep::particle& part){
  /**************************************************************/

  if ( part.charge() != 0 ){
    _chadP[0] += part.px();
    _chadP[1] += part.py();
    _chadP[2] += part.pz();
    _chadP[3] += part.e();
    _nhad[0]++;
  } else {
    _nhadP[0] += part.px();
    _nhadP[1] += part.py();
    _nhadP[2] += part.pz();
    _nhadP[3] += part.e();
    //cout<<"_nhadP[3] ="<<_nhadP[3]<<endl;
  }
  // _hadP[0] += part.px();
  //   _hadP[1] += part.py();
  //   _hadP[2] += part.pz();

  //   _hadE[0] += part.e();
 

}

//*************************************************************************************
void MINDplotter::hadron_direction(fitter& fit) {
  //*************************************************************************************
  
  double normal;
  EVector fitunit;
  for(int i=0; i<2; i++){
    normal = sqrt(pow(_hadP[0],2)+pow(_hadP[1],2)+pow(_hadP[2],2));
  //   if ( _nhits >= 2 ){
    fitunit = fit.get_had_unit()[i];
  
    _haddot[i] = fitunit[0]*(_hadP[0]/normal)
      +fitunit[1]*(_hadP[1]/normal)
      +fitunit[2]*(_hadP[2]/normal);
  
  // } else _haddot = 99;
    
    _hadE[1][i] = fit.get_had_eng()[i];
    _hadE[3][i] = fit.get_had_profile()[i];
    _hrecP[i][0] = fitunit[0]; _hrecP[i][1] = fitunit[1]; _hrecP[i][2] = fitunit[2];
  // else if(_traj->quality("fitted")){
  
    double dxdz = _Th[0][0][_muontraj[i]];
    double dydz = _Th[0][1][_muontraj[i]];
    double mP  = 938.27 * bhep::MeV;
    double mmu = 105.65 * bhep::MeV;
    double mN  = 939.56 * bhep::MeV;
    double Emu = sqrt(pow(1./_qP[0][_muontraj[i]],2) + mmu*mmu); 
    double costhmu = 1./sqrt(1. + dxdz*dxdz + dydz*dydz);
    if(_qP[0][_muontraj[i]] > 0) // assume the QE interaction numubar + p \to n + mu^{+}
      _hadE[2][i] = (mP * Emu + 0.5*(mN*mN - mmu*mmu - mP*mP))/(mP - Emu + costhmu/fabs(_qP[0][_muontraj[i]]));
    else if(_qP[0][_muontraj[i]] < 0) // assume the QE interaction numu + n \to p + mu^{-}
      _hadE[2][i] = (mN * Emu + 0.5*(mP*mP - mmu*mmu - mN*mN))/(mN - Emu + costhmu/fabs(_qP[0][_muontraj[i]]));
    else 
      _hadE[2][i] = 0.;
  }


  ///positions of hits
  Trajectory& traj = fit.get_hadTrajs();

  for (int iHits = 0;iHits < traj.size();iHits++){
	
    const Measurement& meas = traj.node(iHits).measurement();
    
    _XHadPos.push_back(meas.position()[0]);
    _YHadPos.push_back(meas.position()[1]);
    _ZHadPos.push_back(meas.position()[2]);
    
  }

}

//*************************************************************************************
void MINDplotter::max_local_chi2(const Trajectory& traj, const int trajno) {
  //*************************************************************************************

  _m.message("++Finding trajectory local chi max++",bhep::VERBOSE);
  
  size_t nNodes = traj.size();
  double trajMax = 0;
  
  //cout<<traj.quality()<<"  "<<traj.ndof()<<"   "<<nNodes<<"  "<<traj.quality("chi2")<<endl;

  for (size_t iNode = 0;iNode < nNodes;iNode++){

    
    if ( traj.node(iNode).qualitymap().has_key("predicted") )
      trajMax = TMath::Max(trajMax, traj.node(iNode).quality("predicted") );

  }
  
  _Chi[0][trajno] = traj.quality();
  _Chi[1][trajno] = trajMax;
  _Chi[2][trajno] = traj.quality() * _ndf[trajno];
  _m.message("ndof=",_ndf[trajno],"  ",_Chi[0][trajno],"  ",_Chi[2][trajno],bhep::VERBOSE);
}
/*
//****************************************************************************************
void MINDplotter::patternStats1(fitter& Fit, const Trajectory& traj, const int trajno) {
  //****************************************************************************************
  
  _nhits[trajno] = traj.size();
  
  for (int iHits = 0;iHits < _nhits[trajno];iHits++){
    const Measurement& meas = traj.measurement(iHits);

    _XPos[iHits] = meas(iHits)->vector()[0];
    _YPos[iHits] = meas(iHits)->vector()[1];
    _ZPos[iHits] = meas(iHits)->position()[2];

    if (!_patR)
      if ( traj.node(iHits).status("fitted") )
	_hitType[3][trajno]++;
  }

  //Event classifier version.
  //_nhits = Fit.get_nMeas();
  const dict::Key candHit = "inMu";
  const dict::Key engDep = "E_dep";
  int nNode = 0;
  if ( Fit.check_reseed() ) nNode = (int)traj.size()-1;
  bool isMu;
  
  for (int iHits = 0;iHits < _nhits[trajno];iHits++){
  
   
    if (Fit.get_meas(iHits)->name("MotherParticle").compare("mu+")==0
	|| Fit.get_meas(iHits)->name("MotherParticle").compare("mu-")==0){
      isMu = true;
      _mus[trajno][iHits] = true;
      _hitType[0][trajno]++;

            
    }
    else {_mus[trajno][iHits] = false; isMu = false;}
    
    // if ( _fail[trajno] != 7 && Fit.get_meas(iHits)->names().has_key(candHit) ){
    if ( _failEvent != 7 && Fit.get_meas(iHits)->names().has_key(candHit) ){
      if ( Fit.get_meas(iHits)->name(candHit).compare("True")==0 ){//has_key(candHit) ){
	_cand[trajno][iHits] = true;
	_hitType[1][trajno]++;
	_engTraj[trajno] += bhep::double_from_string( Fit.get_meas(iHits)->name(engDep) ) * bhep::GeV;
	if ( isMu ) _hitType[2][trajno]++;
	
	
	if ( traj.node(nNode).status("fitted") && _fail[trajno]!=1 && _fail[trajno]<4 ){	
	  _node[trajno][iHits] = true; _hitType[3][trajno]++; }
	else _node[trajno][iHits] = false;
	if ( Fit.check_reseed() ) nNode--;
	else nNode++;


      }
      else { _node[trajno][iHits] = false; _cand[trajno][iHits] = false; }
    } else if ( _failEvent != 7) { _node[trajno][iHits] = false; _cand[trajno][iHits] = false; }
    //else if ( _fail[trajno] != 7) { _node[trajno][iHits] = false; _cand[trajno][iHits] = false; }

  }
  
  // if (_fail[trajno] != 7){
   if (_failEvent != 7){
    _pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
    _pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
    _pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];

    _engvar[1][trajno] = 0;
    for(int ii=0;ii<_hitType[1][trajno];ii++)
      _engvar[1][trajno] += pow( bhep::double_from_string( Fit.get_traj().node(ii).measurement().name(engDep) ) * bhep::GeV - _engTraj[trajno]/_hitType[1][trajno], 2);
    _engvar[1][trajno] /= (_hitType[1][trajno]-1);
    _engvar[0][trajno] = _engTraj[trajno]/_hitType[1][trajno];
  } else {_engvar[1][trajno] = -1; _engvar[0][trajno] = -1;}
  
  // for (int iclear = _nhits;iclear<300;iclear++){
  //     _mus[iclear] = false; _cand[iclear] = false;
  //     _node[iclear] = false;
  //   }

}
*/


//***********************************************************************/
void MINDplotter::hitBreakUp(fitter& Fit) {
  //***********************************************************************/
  
  ///calculate true muon hits and hadron hits from the event
  _hitTrMu = 0;
  _hitHad = 0;
  const dict::Key hadHit = "inhad";


  /// measurement vector
  std::vector<cluster*>& meas = Fit.GetMeasVec();
  _m.message(" nmeas=   ",Fit.GetNMeas(),bhep::VERBOSE);
  

  ///total measurements/hits in the event
  _hitsInEvent = Fit.GetNMeas();

  ///loop over measurements 
  for(int ih=0; ih<Fit.GetNMeas(); ih++ ){

    //total no of truMu
    // if ( meas[ih]->get_mu_prop() > 0.8 ){
    //  _hitTrMu++;
    if ( meas[ih]->hv("MuonProp").vector()[0] > 0.8 ) { _hitTrMu++;
      _m.message( " TruMu X,Y,Z =",meas[ih]->position()[0],meas[ih]->position()[1],meas[ih]->position()[2],"  No= ",_hitTrMu, bhep::DETAILED);
    }
    
    //total no of hadrons
    if ( meas[ih]->names().has_key(hadHit) ) _hitHad++;

    vector<string> mother_particle = Fit.GetMeas(ih)->get_mother_particle();

    double mu = 0.0;
    cout<<"Mother_particle"<<endl;
    for(int i=0;i<mother_particle.size();i++)
      {
	cout<<mother_particle[i]<<endl;
	_MotherMeas.push_back(mother_particle[i]);

	if(mother_particle[i] == "mu+" || mother_particle[i] == "mu-") mu++; 
	
      }

    _MotherProp.push_back(mu/mother_particle.size());
    


    _XMeas.push_back(Fit.GetMeas(ih)->position()[0]);
    _YMeas.push_back(Fit.GetMeas(ih)->position()[1]);
    _ZMeas.push_back(Fit.GetMeas(ih)->position()[2]);
    _EMeas.push_back(
		     Fit.get_classifier().correctEdep( 
						      Fit.GetMeas(ih)->hv("energy").vector()[0], 
		     Fit.GetMeas(ih)->position()[0], 
		     Fit.GetMeas(ih)->position()[1], 
						      Fit.GetMeas(ih)->position()[2]));
  
    _MuProp.push_back(Fit.GetMeas(ih)->get_mu_prop());

  }
  //// First Angles  
  if ( Fit.GetNMeas() > 6 ){
    //_Theta0 = atan (((_YMeas[0]-_tX[0]) / (_ZMeas[0]-_tX[2])))*180/3.14; 
    _Theta1 = atan (((_YMeas[1]-_YMeas[0]) / (_ZMeas[1]-_ZMeas[0])))*180/3.14;  
    _Theta2 = atan (((_YMeas[3]-_YMeas[2]) / (_ZMeas[3]-_ZMeas[2])))*180/3.14;
    _Theta3 = atan (((_YMeas[5]-_YMeas[4]) / (_ZMeas[5]-_ZMeas[4])))*180/3.14;

    _Delta1 = _Theta2 - _Theta1 ;
    _Delta2 = _Theta3 - _Theta2 ;
  }
  else {
    _Theta1 = 0;
    _Theta2 = 0;
    _Theta3 = 0;
    _Delta1 = 0;
    _Delta2 = 0;
  }

  ////

  _m.message("inside hitBreakUp:: TruMu = ",_hitTrMu,bhep::VERBOSE);
  
  
  ///trajectories from fitter
  std::vector<Trajectory*>& trajs = Fit.get_trajs();
  
  //HitBreakUp for each traj
  for(int tj=0; tj<(int)trajs.size(); tj++){
    
    //
    std::vector<Node*> nhits = trajs[tj]->nodes();
    
    ///loop over hits of each trajectory
    for (int iHits = 0;iHits <(int)trajs[tj]->size(); iHits++){
      if (_fail[tj] != 2){
	
	///candidate hits (all hits are candidates now)
	_hitType[0][tj]++;
	
	///tru muon among cand hits
	if ( nhits[iHits]->measurement().hv("MuonProp").vector()[0] > 0.8) {
	  // if ( nhits[iHits]->get_mu_prop() > 0.8) {
	  _hitType[1][tj]++; 
	  
	  _m.message( "MuInMu X,Y,Z =",nhits[iHits]->measurement().position()[0],nhits[iHits]->measurement().position()[1], nhits[iHits]->measurement().position()[2],"  No= ",_hitType[1][tj], bhep::DETAILED);
	  
	  
	  ///to get true hit index
	  _truHitIndex[0][tj] = iHits;
	  
	}
	///Candidate hit index
	_truHitIndex[1][tj] = iHits;
	
	
	///fitted hits
	// if ( !trajs[tj]->node(iHits).status().has_key("fitted") ) continue;
	if ( trajs[tj]->node(iHits).status("fitted") && _fail[tj]!=1 && _fail[tj]<4 ) _hitType[2][tj]++; 
	// cout<<nhits[iHits]->hv("MuonProp").vector()[0]<<" pos="<<nhits[iHits]->position()[2]<<endl;
      }
    }
    _m.message("traj No= ",tj," ;total InMu=",_hitType[0][tj],"  ;total MuInMu=",_hitType[1][tj], " fitted node=",_hitType[2][tj],bhep::VERBOSE);
  }
  
  
}




/**********************************************************************************/
void MINDplotter::patternStats2(fitter& Fit) {
  /**********************************************************************************/ 

  // position of event hits

  

  ///positions of hits
  std::vector<Trajectory*> &traj = Fit.get_trajs();
  _nallhits=0;

  //loop over trajectories
  for(int i=0; i<(int)traj.size(); i++){
    //if(traj[i]->quality("fitted")){

    if(1){
      _nhits[i] = traj[i]->size();

      //cout<<"nodes="<<traj[i]->size()<<" hits="<<_nhits[i]<<endl;

      /// Sum of hits of trajectory
      _nallhits += traj[i]->size();

      ///temporary vectors
      std::vector<double> vxpos;
      std::vector<double> vypos;
      std::vector<double> vzpos;
      std::vector<double> vedep;
      std::vector<double> vtime;
      std::vector<double> vfitted;
      
      vxpos.clear();
      vypos.clear();
      vzpos.clear();
      vedep.clear();
      vtime.clear();
      vfitted.clear();
      

      for (int iHits = 0;iHits < _nhits[i];iHits++){
	
	const Measurement& meas = traj[i]->node(iHits).measurement();
	
	if(_fail[i] != 2) {
	  vxpos.push_back(meas.position()[0]);
	  vypos.push_back(meas.position()[1]);
	  vzpos.push_back(meas.position()[2]);
	  vedep.push_back(Fit.get_classifier().correctEdep( meas.hv("energy").vector()[0], 
							    meas.vector()[0], meas.vector()[1],
							    meas.vector()[2]));
	  vtime.push_back(meas.hv("hittime").vector()[0]);

	  vfitted.push_back(traj[i]->node(iHits).status("fitted"));
	}
	
	//cout<<"_Xpos ="<<  _XPos[trajno][iHits]<<"  _Ypos ="<<  _YPos[trajno][iHits]<<"   Edep="<<_Edep[trajno][iHits]<<endl;    
	
	//	if (!_patR)
	//if ( traj[i]->node(iHits).status("fitted") )
	//  _hitType[2][i]++;
      }
      
      
      ///Fill the Tree variables for hits
      _XPos.push_back(vxpos);
      _YPos.push_back(vypos);
      _ZPos.push_back(vzpos);
      _Edep.push_back(vedep);
      _HTime.push_back(vtime);
      _NodeFitted.push_back(vfitted);
                  
      ///creat the vector<cluster> from Measurements
      std::vector<double> SelEdep;
      std::vector<Node*> hits = traj[i]->nodes();
      
      //loop over hits
      for (int iHits = 0;iHits < _nhits[i];iHits++){

	if ( _failEvent != 7 && _fail[i] != 2){
	  
	  ///energy of trajectory
	  _engTraj[i] += 
	    Fit.get_classifier().correctEdep(hits[iHits]->measurement().hv("energy").vector()[0]* bhep::MeV,
					     hits[iHits]->measurement().position()[0], hits[iHits]->measurement().position()[1],
					     hits[iHits]->measurement().position()[2]);
	  // cout<<"	_engTraj[i]="<<	_engTraj[i]<<endl;
	  
	    }

      }
      
      
      
      ///pattern Recognition Chi2
      if (_failEvent != 7){
	_pChi[i][0] = Fit.get_classifier().get_PatRec_Chis()[i][0];
	_pChi[i][1] = Fit.get_classifier().get_PatRec_Chis()[i][1];
	_pChi[i][2] = Fit.get_classifier().get_PatRec_Chis()[i][2];
	
	_engvar[1][i] = 0;
	
	
	
	////sum of edep variance from mean edep for all hits of a given track
	for(int ii=0;ii<_nhits[i];ii++)
	  _engvar[1][i] += _fail[i] != 2 ?
	    pow(Fit.get_classifier().correctEdep( (hits[ii]->measurement().hv("energy").vector()[0])*MeV,
						  hits[ii]->measurement().position()[0], hits[ii]->measurement().position()[1], 
						  hits[ii]->measurement().position()[2])
		- _engTraj[i]/_nhits[i], 2) : 1e-6;

	_engvar[1][i] /= (_nhits[i]-1);
	_engvar[0][i] = _engTraj[i]/_nhits[i];
	
	double shortmean = 0.0;
	int Nmean        = 0;
	//if(inMuC.size() > 2){
	//for(int It=0;It<int(inMuC.size());It++){
	if(_nhits[i] > 2){
	  for(int It=0;It<int(_nhits[i]);It++){
	    // first check if position is in the last 2/3 of track
	    // if not increment both iterators and go to the next loop
	   
	    if(It < int(0.3*_nhits[i]) && It > _nhits[i] - 4) continue;

	    // otherwise consider all hits in the detector
	    shortmean += _fail[i] != 2 ? 
	      Fit.get_classifier().correctEdep((hits[It]->measurement().hv("energy").vector()[0])*MeV,
					       hits[It]->measurement().position()[0], hits[It]->measurement().position()[1],
					       hits[It]->measurement().position()[2]) : 0.0;
	    Nmean++;
	    SelEdep.push_back(_fail[i] != 2 ? 
			      Fit.get_classifier().correctEdep((hits[It]->measurement().hv("energy").vector()[0])*MeV,
							       hits[It]->measurement().position()[0], hits[It]->measurement().position()[1],
							       hits[It]->measurement().position()[2]): 1e-6);
	  }
	  
	  sort (SelEdep.begin(), SelEdep.end());
	  double sumlow=0.0, sumhigh=0.0;
	  int N = int(SelEdep.size()/2);
	  for (int jj=0; jj<(int)SelEdep.size(); jj++){
	    if(jj < N)
	      sumlow += SelEdep[jj];
	    else
	      sumhigh += SelEdep[jj];
	  }
	  _engvar[2][i] = N > 0 ? sumlow/double(N) : 0.0;
	  _engvar[3][i] = (int)SelEdep.size() >= N ? sumhigh/double(SelEdep.size() - N) : 1.0;
	} else {
	  _engvar[2][i] = 0.0; _engvar[3][i] = 1.0;
	}
	
	///clear the vectors
	hits.clear();
	SelEdep.clear();
	
      } else {_engvar[1][i] = -1; _engvar[0][i] = -1; 
	_engvar[2][i] = 0; _engvar[3][i] = 1;}
      
      /// inMuC.clear();
    }
  }
}




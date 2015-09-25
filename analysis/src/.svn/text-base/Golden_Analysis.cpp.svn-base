// Implementation of class Golden_Analysis.

#include "Golden_Analysis.h"
#include <TRandom3.h>
//#include <CLHEP/Random/RandGauss.h>
//#include <CLHEP/Random/RandFlat.h>

/********************************************/
Golden_Analysis::Golden_Analysis(const bhep::gstore& store)
/********************************************/
{
  //Initialise all cuts and relevant histograms.
  _filebase = store.fetch_sstore("filebase");

  std::string outName = store.fetch_sstore("outFile");
  _OutFile = new TFile( outName.c_str(), "recreate" );

  _intType = store.fetch_sstore("intType");
  _anType = store.fetch_istore("anType");
  _beamCharge = store.fetch_istore("bChar");
  _doSmear = store.fetch_istore("doSmear");

  define_histos( store );

  define_cut_levels( store );

  _firstF = store.fetch_istore("firstF");
  _lastF = store.fetch_istore("lastF");

  if ( _doSmear ){
    long seed = (long)store.fetch_dstore("seed");
    _ranGen = TRandom3( seed );
    _smearRoot = store.fetch_dstore("rootP");
    _smearConst = store.fetch_dstore("constP");
    _angA = store.fetch_dstore("angA");
    _angB = store.fetch_dstore("angB");
  }

  _doVary = store.fetch_istore("doVary");

  if ( _doVary > 0 )
    define_variables( store );

  _detX = store.fetch_dstore("detX");
  _detY = store.fetch_dstore("detY");
  _WLSAtten = store.fetch_dstore("WLSatten");

  isMeson = true;
}

/********************************************/
void Golden_Analysis::execute()
/********************************************/
{
  //Do the analysis!!
  std::string fileName;
  bool includeEvent;
  
  for (int ifile = _firstF;ifile <= _lastF;ifile++){
    
    fileName = _filebase+bhep::to_string( ifile )+".root";
    std::cout<<fileName<<std::endl;
    TFile *f1 = new TFile( fileName.c_str(), "read" );
    if(f1->IsZombie()) continue;
    TTree *data;
    f1->GetObject("tree", data);

    set_branches( *data );

    for (int ient = 0;ient < data->GetEntries();ient++){
      
      data->GetEntry( ient );
      if(ient%1000 == 0) std::cout<<"Event = "<<ient<<std::endl;
      if ( _doVary > 0 )
	includeEvent = include_event();
      else includeEvent = true;
      isMeson = true;
      if(_npi[0] + _npi[1] + _npi[2] + _nk[0] + _nk[1] + _nk[2] == 0)
	isMeson = false;
      //avoids inverse muon decays which can get into NC sample
      //when generating with Nuance.
      if ( _intType == "NC" && _charge[0] != 0 )
	std::cout << "Avoiding CC in NC sample" <<std::endl;
      else {
	if ( includeEvent ){

	  _alltracks1->Fill( _enu/1000 );
	  IntHist[_truInt]->Fill( _enu/1000 );
	  _ifail->Fill(_fail);
	  fill_tracks2();
	  if ( _intType == "nu_mu" ) _ally->Fill( 1 - sqrt( pow(1.0/_mom[0],2)+pow(105.66,2) )/_enu );
	  
	  if ( _suc )
	    analyse_success();
	}
      }
    }
    
    delete data;
    delete f1;

  }
  
  // if ( _intType == "nu_mu" ){
    _Eff->Divide( _alltracks1 );
    _recEff->Divide( _alltracks2 );
    _Ey->Divide( _ally );
    // }
  _back->Divide( _alltracks1 );
  _recBack->Divide( _alltracks2 );
  for(int i=0;i<10;i++){
    IntHist[i]->Divide( _alltracks1 );
    IntEff[i]->Divide( IntHist[i] );
    Intback[i]->Divide( IntHist[i] );
  }
  _EffNM->Divide( _alltracksNM1 );
  _backNM->Divide( _alltracksNM1 );
  _By->Divide( _ally );
  _rcuts->Scale(1./_alltracks1->Integral());
  _pcuts->Scale(1./_alltracks1->Integral());
  _rcuts_back->Scale(1./_alltracks1->Integral());
  _pcuts_back->Scale(1./_alltracks1->Integral());
  _rcuts_sig->Scale(1./_alltracks1->Integral());
  _pcuts_sig->Scale(1./_alltracks1->Integral());
  _LCC4_sig->Scale(1./_LCC4_sig->Integral());
  _LCC4_back->Scale(1./_LCC4_back->Integral());
  _ifail->Scale(1./_alltracks1->Integral());
}

/********************************************/
void Golden_Analysis::finalize()
/********************************************/
{
  //Finalise outputs
  _OutFile->Write();
  _OutFile->Close();
  delete _OutFile;
  
  _likeCC->Close();
  _likeNC->Close();
  delete _likeCC;
  delete _likeNC;
  
  if ( _doVary > 0 ){
    _xsecF->Close();
    delete _xsecF;
  }

}

/********************************************/
void Golden_Analysis::fill_tracks2()
/********************************************/
{
  //Fill a column of the 2D matrix for normalisation.
  double incr = 0.0;
  int bin;
  while ( incr < _alltracks2->GetXaxis()->GetXmax() )
    {
      _alltracks2->Fill( incr, _enu/1000 );

      bin = _alltracks2->GetXaxis()->FindBin( incr );
      incr += _alltracks2->GetXaxis()->GetBinWidth( bin );
    }

}

/********************************************/
void Golden_Analysis::define_histos(const bhep::gstore& store)
/********************************************/
{
  //Define the required output histograms.
  _maxE = store.fetch_dstore("truEmax");
  int ntBins = store.fetch_istore("truBins");
  double trEdge[ntBins+1], recEdge[ntBins+2];
  define_rec_binning( store, recEdge, trEdge );

  //if ( _intType == "nu_mu" ){

  _Eff = new TH1F("Eff",";True Neutrino Energy;Fractional Efficiency",(int)_maxE,0,_maxE);
  _Eff->Sumw2();

  _EffNM = new TH1F("EffNM",";True Neutrino Energy;Fractional Efficiency",(int)_maxE,0,_maxE);
  _EffNM->Sumw2();
  
  _Ey = new TH1F("Ey","Efficiency as function of inelasticity",100,0,1);
  _Ey->Sumw2();
  
  _By = new TH1F("By","Background as function of inelasticity",100,0,1);
  _By->Sumw2();
  
  _recEff = new TH2F("recEff","Signal response",ntBins+1, recEdge, ntBins, trEdge);
  _recEff->Sumw2();
  
  _recmunu = new TH2F("recmunu",";Rec. #nu Energy (GeV/c); Rec. #mu Energy (GeV/c)",
		      (int)_maxE, 0, _maxE, (int)_maxE, 0, _maxE);
  _recmunu->Sumw2();
  
  _recQtnu = new TH2F("recQtnu",";Rec. #nu Energy (GeV/c); Q_{t} (GeV/c)",
		      (int)_maxE, 0, _maxE, (int)_maxE, 0, _maxE/10.0);
  _recQtnu->Sumw2();
  

  _dir_back = new TH2D("dir_back",";p_{x}/p_{z};p_{y}/p_{z}", 100, -2.5, 2.5, 100, -2.5, 2.5);
  _dir_back->Sumw2();
  _dir_sig = new TH2D("dir_sig",";p_{x}/p_{z};p_{y}/p_{z}", 100, -2.5, 2.5, 100, -2.5, 2.5);
  _dir_sig->Sumw2();
  _dirMC_back = new TH2D("dirMC_back",";p_{x}/p_{z};p_{y}/p_{z}", 100, -2.5, 2.5, 100, -2.5, 2.5);
  _dirMC_back->Sumw2();
  _dirMC_sig = new TH2D("dirMC_sig",";p_{x}/p_{z};p_{y}/p_{z}", 100, -2.5, 2.5, 100, -2.5, 2.5);
  _dirMC_sig->Sumw2();

  _ally = new TH1F("tracksy","All events inelasticity",100,0,1);
  _ally->Sumw2();
  
  _pmu = new TH1F("p_mu",";Reconstructed Muon Momentum (GeV/c)",(int(_maxE*5)),0,_maxE);
  _pmu->Sumw2();

  _pmutfail = new TH1F("p_mutfail",";True Muon Momentum (GeV/c)",(int(_maxE*5)),0,_maxE);
  _pmutfail->Sumw2();

  _pmufail = new TH1F("p_mufail",";Reconstructed Muon Momentum (GeV/c)",(int(_maxE*5)),-_maxE,_maxE);
  _pmufail->Sumw2();
  
  /*} else {

    _Eff = NULL;
    _recEff = NULL;
    _recmunu = NULL;
    _recQtnu = NULL;
    _Ey = NULL;
    _By = NULL;
    _ally = NULL;
    _pmu = NULL;

    }*/
  _rcuts = new TH1F("cuts_reject", "Events rejected by cuts",8, 0.0, 8.0);
  _rcuts->Sumw2();
  _pcuts = new TH1F("cuts_passed", "Events passed by cuts",9, 0.0, 9.0);
  _pcuts->Sumw2();

  _rcuts_sig = new TH1F("cuts_reject_sig", "Events rejected by cuts",8, 0.0, 8.0);
  _rcuts_sig->Sumw2();
  _pcuts_sig = new TH1F("cuts_passed_sig", "Events passed by cuts",9, 0.0, 9.0);
  _pcuts_sig->Sumw2();

  _rcuts_back = new TH1F("cuts_reject_back", "Events rejected by cuts",8, 0.0, 8.0);
  _rcuts_back->Sumw2();
  _pcuts_back = new TH1F("cuts_passed_back", "Events passed by cuts",9, 0.0, 9.0);
  _pcuts_back->Sumw2();

  _zFid   = new TH1F("zFid",";Vertex Z-position (m)", 140, -70, 70);
  _zFidcut = new TH1F("zFidcut",";Vertex Z-position (m)", 120, -60, 60);
  _Lqp    = new TH1F("Lqp", ";L_{q/p};Occupancy",10,-5, 5);
  _Lqp_back = new TH1F("Lqp_back", ";L_{q/p};Occupancy",10,-5, 5);
  _LCC1   = new TH1F("LCC1", ";L_1;Occupancy",35, -15, 20);
  _LCC2   = new TH1F("LCC2", ";L_2;Occupancy",35, -15, 20);  
  _LCC3   = new TH1F("LCC3", ";L_3;Occupancy",35, -15, 20);
  _LCC4_sig    = new TH1F("LCC4_sig", ";L_{4};Occupancy",35, -15, 20);
  _LCC4_back   = new TH1F("LCC4_back", ";L_{4};Occupancy",35, -15, 20);
  _ifail = new TH1F("ifail",";Fail Code; Fraction of Events",8,-0.5,7.5);
  _fitprop = new TH1F("fitprop",";Clusters in Fit/Clusters in Track",120,0.0,1.2);

  _back = new TH1F("back","true eng background",(int)_maxE,0,_maxE);
  _back->Sumw2();
  _back_pion = new TH2F("back_pion","true eng background",(int)_maxE,0,_maxE,3,-1.0,2.0);
  _back_pion->Sumw2();
  _back_kaon = new TH2F("back_kaon","true eng background",(int)_maxE,0,_maxE,3,-1.0,2.0);
  _back_kaon->Sumw2();

  _backNM = new TH1F("backNM","true eng background",(int)_maxE,0,_maxE);
  _backNM->Sumw2();

  _recBack = new TH2F("recBack","Background response",ntBins+1, recEdge, ntBins, trEdge);
  _recBack->Sumw2();

  _alltracks1 = new TH1F("tracks","true eng background",(int)_maxE,0,_maxE);
  _alltracks1->Sumw2();

  for(int i=0; i<10; i++){
    TString tempname="nutracks_";
    IntHist.push_back(new TH1F(tempname + (long int)i,";True Neutrino Energy;Fractional Content",
			       2*(int)_maxE,0,_maxE));
    IntHist.back()->Sumw2();
    tempname="nuEff_";
    IntEff.push_back(new TH1F(tempname + (long int)i,";True Neutrino Energy;Fractional Content",
			       2*(int)_maxE,0,_maxE));
    IntEff.back()->Sumw2();
    tempname="nuback_";
    Intback.push_back(new TH1F(tempname + (long int)i,";True Neutrino Energy;Fractional Content",
			       2*(int)_maxE,0,_maxE));
    Intback.back()->Sumw2();
  }

  _alltracksNM1 = new TH1F("tracksNM",";True Neutrino Energy",(int)_maxE,0,_maxE);
  _alltracksNM1->Sumw2();

  _alltracks2 = new TH2F("tracks2D","All event response",ntBins+1, recEdge, ntBins, trEdge);
  _alltracks2->Sumw2();

  _qCharge = new TH1F("qCharge",";#delta c/c",100, -2, 1);
  _qCharge_back = new TH1F("qCharge_back",";#delta c/c",100, -2, 1);

  //Likelihoods.
  std::string lCCname = store.fetch_sstore("cclike");
  std::string lNCname = store.fetch_sstore("nclike");

  _likeCC = new TFile( lCCname.c_str(), "read" );
  _likeNC = new TFile( lNCname.c_str(), "read" );

  _likeCC->GetObject("trHits", _hitCC); _likeNC->GetObject("trHits", _hitNC);
  _likeCC->GetObject("pID", _errSig); _likeNC->GetObject("pnoID", _errBack);
  //err_nc should be for background so combine and renormalise.
  _errBack->Add( (TH1F*)_likeCC->Get("pnoID") );
  _errBack->Scale( 1.0/_errBack->Integral() );

  if ( _anType == 1 ){
    _fracCC = NULL; _fracNC = NULL;
    _meanCC = NULL; _meanNC = NULL;
  } else {
    _likeCC->GetObject("engfracC", _fracCC); _likeNC->GetObject("engfracC", _fracNC);
    _likeCC->GetObject("trMeanC2", _meanCC); _likeNC->GetObject("trMeanC2", _meanNC);
  }
  //
}

/********************************************/
void Golden_Analysis::define_rec_binning(const bhep::gstore& store, double* recB, double* truB)
/********************************************/
{
  //Sort out the bin widths for response matrices.
  vector<double> edges = store.fetch_vstore("bEdges");

  vector<double>::iterator bIt;
  int counter = 0;
  for (bIt = edges.begin();bIt < edges.end();bIt++, counter++)
    {
      truB[counter] = (*bIt);
      recB[counter] = (*bIt);
      if ( bIt == edges.end()-1 ) recB[counter+1] = (*bIt)+1;
    }
  
}

/********************************************/
void Golden_Analysis::define_variables(const bhep::gstore& store)
/********************************************/
{
  //Get the information required for xsec study.
  _Vtype = store.fetch_istore("vtype");

  if ( !_doSmear ){
    long seed = (long)store.fetch_dstore("seed");
    _ranGen = TRandom3( seed );
  }

  std::string xsecFile = store.fetch_sstore("xsecFile");
  _xsecF = new TFile( xsecFile.c_str(), "read" );

  _xsecF->GetObject("xsec", _xsecErr);
}

/********************************************/
bool Golden_Analysis::include_event()
/********************************************/
{
  //Check if this event should be included in the sample.
  if ( _doVary == 1 && _truInt <= _Vtype ) return true;
  if ( _doVary == 2 && _truInt >  _Vtype ) return true;
  
  int bin = _xsecErr->FindBin( _enu/1000 );
  double xsErr = _xsecErr->GetBinContent( bin );
  
  if ( _doVary == 1 ){

    if ( _ranGen.Rndm() <= 1-xsErr ) return true;
    else return false;

  } else if ( _doVary == 2 ){

    if ( _ranGen.Rndm() <= 1/(xsErr+1) ) return true;
    else return false;

  }

  return true;
}

/********************************************/
void Golden_Analysis::define_cut_levels(const bhep::gstore& store)
/********************************************/
{
  //set the cut levels required.
  _edges = store.fetch_vstore("edges");

  _logQPmin = store.fetch_dstore("minQP");

  _maxOver = store.fetch_dstore("maxOver");
  
  _minNodes = store.fetch_dstore("nodeprop");
  
  _likeCuts = store.fetch_vstore("NClike");

  _kMinE = store.fetch_dstore("kinE");
  
  _minQT = store.fetch_dstore("QTmin");
  
  _momGrad = store.fetch_dstore("pGrad");
  
  _quadFitCuts = store.fetch_vstore("quadCut");
  
  _translongCuts = store.fetch_vstore("tranlong");

}

/********************************************/
void Golden_Analysis::set_branches(TTree& inTree)
/********************************************/
{
  //Set the relevant input branches.
  inTree.SetBranchStatus("*", 0);
  inTree.SetBranchStatus("Fitted", 1); inTree.SetBranchAddress("Fitted", &_suc);
  inTree.SetBranchStatus("NeuEng", 1); inTree.SetBranchAddress("NeuEng", &_enu);
  inTree.SetBranchStatus("Charge", 1); inTree.SetBranchAddress("Charge", &_charge);
  inTree.SetBranchStatus("TrueInteraction", 1); inTree.SetBranchAddress("TrueInteraction", &_truInt);
  inTree.SetBranchStatus("Fail",1); inTree.SetBranchAddress("Fail",&_fail);
  inTree.SetBranchStatus("Momentum", 1); inTree.SetBranchAddress("Momentum", &_mom);
  inTree.SetBranchStatus("HitBreakDown", 1); inTree.SetBranchAddress("HitBreakDown", &_trHits);
  inTree.SetBranchStatus("EngDeposit", 1); inTree.SetBranchAddress("EngDeposit", &_depEng);
  inTree.SetBranchStatus("hadEng", 1); inTree.SetBranchAddress("hadEng", &_hadEng);
  inTree.SetBranchStatus("hadronMom", 1); inTree.SetBranchAddress("hadronMom", &_hadMom);
  inTree.SetBranchStatus("FitChiInfo", 1); inTree.SetBranchAddress("FitChiInfo", &_maxLoc);
  inTree.SetBranchStatus("Direction", 1); inTree.SetBranchAddress("Direction", &_dir);
  inTree.SetBranchStatus("CandHits", 1); inTree.SetBranchAddress("CandHits", &_inCand);
  inTree.SetBranchStatus("HadronNodes", 1); inTree.SetBranchAddress("HadronNodes", &_hadH);
  inTree.SetBranchStatus("NoHits", 1); inTree.SetBranchAddress("NoHits", &_nhits);
  inTree.SetBranchStatus("XPositions", 1); inTree.SetBranchAddress("XPositions", &_xPos);
  inTree.SetBranchStatus("YPositions", 1); inTree.SetBranchAddress("YPositions", &_yPos);
  inTree.SetBranchStatus("ZPositions", 1); inTree.SetBranchAddress("ZPositions", &_zPos);
  inTree.SetBranchStatus("npi",1); inTree.SetBranchAddress("npi", &_npi);
  inTree.SetBranchStatus("nk" ,1); inTree.SetBranchAddress("nk",  &_nk );

}

/********************************************/
void Golden_Analysis::analyse_success()
/********************************************/
{
  //Fill graphs if accepted.
  if ( _charge[1] == _beamCharge && _intType == "nu_mu" ){
    _pcuts->Fill("Reconstruction Success",1);
    _pcuts_sig->Fill("Reconstruction Success",1);
    if ( fiducial( 1 ) )
      if ( full_cut( 1 ) ){
	_Eff->Fill( _enu/1000 );
	IntEff[_truInt]->Fill( _enu/1000 );
	if(!isMeson){
	  _alltracksNM1->Fill(_enu/1000);
	  _EffNM->Fill( _enu/1000 );
	}
	if ( _recEnu > _maxE*1000 ) _recEnu = _maxE*1000;//Put into overflow bin.
	_recEff->Fill( _recEnu/1000, _enu/1000 );
	_Ey->Fill( 1 - sqrt( pow(1.0/_mom[0],2)+pow(105.66,2) )/_enu );
	_pmu->Fill(1.0/_mom[1]/1000);
	_dir_sig->Fill(_dir[1][0], _dir[1][1]);
	_dirMC_sig->Fill(_dir[0][0], _dir[0][1]);
      }
    
  } else if ( _charge[1] == -_beamCharge){
    _pcuts->Fill("Reconstruction Success",1);
    _pcuts_back->Fill("Reconstruction Success",1);
    if ( fiducial(-1 ) )
      if ( full_cut( -1 ) ){
	_back->Fill( _enu/1000 );
	Intback[_truInt]->Fill( _enu/1000 );
	if(!isMeson) {
	  _backNM->Fill( _enu/1000 );
	  _alltracksNM1->Fill(_enu/1000);
	}
	if ( _recEnu > _maxE*1000 ) _recEnu = _maxE*1000;//Put into overflow bin.
	_recBack->Fill( _recEnu/1000, _enu/1000 );
	if ( _intType == "nu_mu" ){
	  _By->Fill( 1 - sqrt( pow(1.0/_mom[0],2)+pow(105.66,2) )/_enu );
	  _pmutfail->Fill(1./_mom[0]/1000);
	  _pmufail->Fill(1./_mom[1]/1000);
	}
	_dir_back->Fill(_dir[1][0], _dir[1][1]);
	_dirMC_back->Fill(_dir[0][0], _dir[0][1]);
	if ( _npi[0] > 0) _back_pion->Fill( _enu/1000, +1);
	if ( _npi[1] > 0) _back_pion->Fill( _enu/1000,  0);
	if ( _npi[2] > 0) _back_pion->Fill( _enu/1000, -1);
	if ( _nk[0]  > 0) _back_kaon->Fill( _enu/1000, +1);
	if ( _nk[1]  > 0) _back_kaon->Fill( _enu/1000,  0);
	if ( _nk[2]  > 0) _back_kaon->Fill( _enu/1000, -1);
      }   
  }
  
}

/********************************************/
bool Golden_Analysis::fiducial(int ID)
/********************************************/
{
  //controls the edge rejection.
  int i1 = -1, i2 = _nhits;
  double x1, xend, y1, yend, z1, zend;

  do {
    i1++;

    if ( _inCand[i1] ){
      x1 = _xPos[i1];
      y1 = _yPos[i1];
      z1 = _zPos[i1];
    }

  } while ( i1 < _nhits && !_inCand[i1] );
  _zFid->Fill(z1/1000);
  if ( z1 > _edges[0] ) {
    _zFidcut->Fill(z1/1000);
    _rcuts->Fill("Fiducial",1); 
    if(ID==1) _rcuts_sig->Fill("Fiducial",1);
    else _rcuts_back->Fill("Fiducial",1);
    return false; 
  }
  _pcuts->Fill("Fiducial",1);
  if(ID==1) _pcuts_sig->Fill("Fiducial",1);
  else _pcuts_back->Fill("Fiducial",1);
  
  
//   do {
//     i2--;

//     if ( _inCand[i2] ){
//       xend = _xPos[i2];
//       yend = _yPos[i2];
//       zend = _zPos[i2];
//     }

//   } while ( i2 >= _nhits && !_inCand[i2] );

//   //Accepted in fiducial;
//   if ( zend > _edges[1] && z1 > _edges[0]) return false;
//   if ( fabs(xend) > _edges[3] && fabs(x1) > _edges[2] )
//     if ( xend/fabs(xend) == x1/fabs(x1) ) return false;
//   if ( fabs(yend) > _edges[5] && fabs(y1) > _edges[4] )
//     if ( yend/fabs(yend) == y1/fabs(y1) ) return false;

  return true;
}

/********************************************/
bool Golden_Analysis::full_cut(int ID)
/********************************************/
{
  //Sequencially perform cuts.
  if ( fabs( _mom[2]/_mom[1]) >= _errSig->GetXaxis()->GetXmax() ) {
      return false;}
  else {
    
    double err_cc = _errSig->GetBinContent( _errSig->FindBin( fabs( _mom[2]/_mom[1] ) ) );
    double err_nc = _errBack->GetBinContent( _errBack->FindBin( fabs( _mom[2]/_mom[1] ) ) );
    double errPar = log( err_cc/err_nc );
    
    ID==1? _Lqp->Fill(errPar): _Lqp_back->Fill(errPar);
    
    if ( fabs(1.0/(1000*_mom[1])) > _maxOver*_maxE ){// 
      _rcuts->Fill("Max Momentum",1); 
      if(ID==1) _rcuts_sig->Fill("Max Momentum",1); 
      else _rcuts_back->Fill("Max Momentum",1); 
      return false;
    }
    else {
      _pcuts->Fill("Max Momentum",1);
      if(ID==1) _pcuts_sig->Fill("Max Momentum",1); 
      else _pcuts_back->Fill("Max Momentum",1); 
    }
    _fitprop->Fill((double)_trHits[3]/(double)_trHits[1]);
    if ( (double)_trHits[3]/(double)_trHits[1] < _minNodes ){ 
      _rcuts->Fill("Fitted proportion",1); 
      if(ID==1) _rcuts_sig->Fill("Fitted proportion",1);
      else _rcuts_back->Fill("Fitted proportion",1);
      return false; 
    } // 
    else{
      _pcuts->Fill("Fitted proportion",1);
      if(ID==1) _pcuts_sig->Fill("Fitted proportion",1);
      else _pcuts_back->Fill("Fitted proportion",1);
    }
    if ( errPar <= _logQPmin ) { // 
      _rcuts->Fill("Track quality",1); 
      if(ID==1) _rcuts_sig->Fill("Track quality",1); 
      else _rcuts_back->Fill("Track quality",1); 
      return false;}
    else {
      _pcuts->Fill("Track quality",1);
      if(ID==1) _pcuts_sig->Fill("Track quality",1); 
      else _pcuts_back->Fill("Track quality",1); 
    }

    if ( !quad_and_transLong(ID) )
      return false; 
    
    if ( !NC_likelihood(ID) ){ // 
      _rcuts->Fill("CC Selection",1); 
      if(ID==1) _rcuts_sig->Fill("CC Selection",1); 
      else _rcuts_back->Fill("CC Selection",1); 
      return false;}
    else{
      _pcuts->Fill("CC Selection",1);
      if(ID==1) _pcuts_sig->Fill("CC Selection",1); 
      else _pcuts_back->Fill("CC Selection",1); 
    }

    if ( !kin_cut( ID ) ) { // 
      _rcuts->Fill("Kinematic",1); 
      if(ID==1) _rcuts_sig->Fill("Kinematic",1); 
      else _rcuts_back->Fill("Kinematic",1); 
      return false;}
    else{ _pcuts->Fill("Kinematic",1);
      if(ID==1) _pcuts_sig->Fill("Kinematic",1); 
      else _pcuts_back->Fill("Kinematic",1); 
    }
  }
  
  return true;
}

/********************************************/
bool Golden_Analysis::NC_likelihood(int ID)
/********************************************/
{
  //Neutral current likelihood rejection.
  double CChit, NChit;
  double CCfrac, NCfrac;
   double CCmean, NCmean;
  double logPar;

  if ( _anType == 1 ){

    calc_hit_proportion();
    
    if ( _trHits[1] >= _hitCC->GetXaxis()->GetXmax() )
      return true;
    else {

      logPar = calc_hitlog();
      
      if ( logPar > _likeCuts[0] ) return true;
    }

  } else {

    if ( _trHits[1] >= _hitCC->GetXaxis()->GetXmax() ){
      calc_hit_proportion();
      return true;
    } else {    

      calc_visible_eng();

      if ( _trHits[1] < _hitCC->GetXaxis()->GetXmax() &&
	   _trajEng/(double)_trHits[1] >= _meanCC->GetXaxis()->GetXmax() && _trajEng/_visEng > 0.999 ){
	
	logPar = calc_hitlog();
	_LCC1->Fill(logPar);
	if ( logPar > _likeCuts[0] ) return true;

      } else if ( _trHits[1] < _hitCC->GetXaxis()->GetXmax() &&
		  _trajEng/(double)_trHits[1] >= _meanCC->GetXaxis()->GetXmax() && _trajEng/_visEng <= 0.999 ){
	
	logPar = calc_hitfraclog();
	_LCC2->Fill(logPar);
	if ( logPar > _likeCuts[1] ) return true;
	
      } else if ( _trHits[1] < _hitCC->GetXaxis()->GetXmax() &&
		  _trajEng/(double)_trHits[1] < _meanCC->GetXaxis()->GetXmax() && _trajEng/_visEng > 0.999 ){
	
	logPar = calc_hitmeanlog();
	_LCC3->Fill(logPar);
	if ( logPar > _likeCuts[2] ) return true;
	
      } else if ( _trHits[1] < _hitCC->GetXaxis()->GetXmax() &&
		  _trajEng/(double)_trHits[1] < _meanCC->GetXaxis()->GetXmax() && _trajEng/_visEng <= 0.999 ){
	
	logPar = calc_all3log();
	ID==1? _LCC4_sig->Fill(logPar):_LCC4_back->Fill(logPar);
	if ( logPar > _likeCuts[3] ) return true;
      }

    }

  }
  
  
  return false;
}

/********************************************/
double Golden_Analysis::calc_hitlog()
/********************************************/
{
  double CChit, NChit;
  CChit = _hitCC->GetBinContent( _hitCC->FindBin( _trHits[1] ) );
  NChit = _hitNC->GetBinContent( _hitNC->FindBin( _trHits[1] ) );

  return log( CChit/NChit );
}

/********************************************/
void Golden_Analysis::calc_hit_proportion()
/********************************************/
{
  //Calculate correct proportion of candidate/event when no eDep.
  _visEng = 0.1;
  _trajEng = 0.01;//avoids pole/incorrect quasi identification.

  int nChit = 0, nEhit = 0;
  for (int Ej=0;Ej<_nhits;Ej++){
	    
    if ( _inCand[Ej] || _hadH[Ej] ) nEhit++;
    if ( _inCand[Ej] ) nChit++;
    
  }
  _trHitProp = (double)nChit/(double)nEhit;

}

/********************************************/
void Golden_Analysis::calc_visible_eng()
/********************************************/
{
  //Do energy correction and calculate trajectory and event visible energy.
  _visEng = 0;
  _trajEng = 0;
  double corrEng;
  int nChit = 0, nEhit = 0;
  for (int Ej=0;Ej<_nhits;Ej++){
	    
    corrEng = correctEdep(_depEng[Ej], _xPos[Ej], _yPos[Ej]);
	    
    if ( _inCand[Ej] || _hadH[Ej] ) {
      _visEng += corrEng;
      nEhit++;
    }
    if ( _inCand[Ej] ) {
      _trajEng += corrEng;
      nChit++;
    }
  }
  _trHitProp = (double)nChit/(double)nEhit;

}

/********************************************/
double Golden_Analysis::correctEdep(double edep, double X, double Y)
/********************************************/
{
  double corrEng;
  
  double sum1 = exp( -(_detX/2-fabs(X))/_WLSAtten ) + exp( -(_detX/2+fabs(X))/_WLSAtten )
    + exp( -(_detY/2-fabs(Y))/_WLSAtten ) + exp( -(_detY/2+fabs(Y))/_WLSAtten );

  corrEng = 4*edep/sum1;
  
  return corrEng;
}

/********************************************/
double Golden_Analysis::calc_hitfraclog()
/********************************************/
{
  double CChit, NChit;
  double CCfrac, NCfrac;
  CChit = _hitCC->GetBinContent( _hitCC->FindBin( _trHits[1] ) );
  NChit = _hitNC->GetBinContent( _hitNC->FindBin( _trHits[1] ) );
  CCfrac = _fracCC->GetBinContent( _fracCC->FindBin( _trajEng/_visEng ) );
  NCfrac = _fracNC->GetBinContent( _fracNC->FindBin( _trajEng/_visEng ) );

  return log( (CChit*CCfrac)/(NChit*NCfrac) );
}

/********************************************/
double Golden_Analysis::calc_hitmeanlog()
/********************************************/
{
  double CChit, NChit;
  double CCmean, NCmean;
  CChit = _hitCC->GetBinContent( _hitCC->FindBin( _trHits[1] ) );
  NChit = _hitNC->GetBinContent( _hitNC->FindBin( _trHits[1] ) );
  CCmean = _meanCC->GetBinContent( _meanCC->FindBin( _trajEng/(double)_trHits[1] ) );
  NCmean = _meanNC->GetBinContent( _meanNC->FindBin( _trajEng/(double)_trHits[1] ) );

  return log( (CChit*CCmean)/(NChit*NCmean) );
}

/********************************************/
double Golden_Analysis::calc_all3log()
/********************************************/
{
  double CChit, NChit;
  double CCfrac, NCfrac;
  double CCmean, NCmean;
  CChit = _hitCC->GetBinContent( _hitCC->FindBin( _trHits[1] ) );
  NChit = _hitNC->GetBinContent( _hitNC->FindBin( _trHits[1] ) );
  CCfrac = _fracCC->GetBinContent( _fracCC->FindBin( _trajEng/_visEng ) );
  NCfrac = _fracNC->GetBinContent( _fracNC->FindBin( _trajEng/_visEng ) );
  CCmean = _meanCC->GetBinContent( _meanCC->FindBin( _trajEng/(double)_trHits[1] ) );
  NCmean = _meanNC->GetBinContent( _meanNC->FindBin( _trajEng/(double)_trHits[1] ) );

  return log( (CChit*CCmean*CCfrac)/(NChit*NCmean*NCfrac) );
}

/********************************************/
bool Golden_Analysis::kin_cut(int ID)
/********************************************/
{
  //Reconstruct energy of neutrino.
  //then use it to do kinematical cuts.
  double Esig, QT, dTh;
  if ( _trHitProp == 1 ){
    //Want quasi eng rec.
    if ( ID == 1 && _beamCharge == 1 )
      _recEnu = quasi_nubar_prot();
    else if ( ID == 1 && _beamCharge == -1 )
      _recEnu = quasi_nu_neut();
    else if ( ID == -1 && _beamCharge == 1 )
      _recEnu = quasi_nu_neut();
    else if ( ID == -1 && _beamCharge == -1 )
      _recEnu = quasi_nubar_prot();
    _recmunu->Fill(_recEnu/1000,fabs(0.001/_mom[1]));
    if ( _recEnu/1000 <= _kMinE ) return true;
    else if ( fabs(1.0/_mom[1]) >= _momGrad*_recEnu ) return true;
  } else {
    //Standard reconstruction. Smear or hadRec.
    if ( _doSmear ){

      if ( ID == -1 && _intType == "NC" )
	_hadEng[0] -= sqrt( pow(1.0/_mom[1],2)+pow(139.57,2) );

      Esig = sqrt( pow( _smearRoot*sqrt( _hadEng[0]/1000 ), 2 )
		   + pow(  _smearConst*_hadEng[0]/1000, 2 ) );

      _recEnu = sqrt( pow(1.0/(_mom[1]*1000),2)+pow(0.10566,2) )
	+ _hadEng[0]/1000 + _ranGen.Gaus( 0, Esig );
      _recEnu *= 1000;

      dTh = sqrt( pow(_angA/sqrt(_hadEng[0]/1000),2) + pow(_angB/(_hadEng[0]/1000),2) );
      QT = calcQt( _ranGen.Gaus( 0, dTh ) );

      _recmunu->Fill(_recEnu/1000,fabs(0.001/_mom[1]));
      _recQtnu->Fill(_recEnu/1000,QT);
      if ( _recEnu/1000 > _kMinE-2.0 && QT <= _minQT ) return false;
      else if ( fabs(1.0/_mom[1]) >= _momGrad*_recEnu || _recEnu/1000 <= _kMinE )
	return true;
      
    } else {
      std::cout << "You twat! You haven't written the rec code yet!"<< std::endl
		<< "Use the smear!" << std::endl;
      exit(0);
    }
  }

  return false;
}

/********************************************/
double Golden_Analysis::quasi_nubar_prot()
/********************************************/
{
  //Quasi formula: nubar p scattering.
  double Eng;
  Eng = 938.27*sqrt( pow(1.0/_mom[1],2)+pow(105.66,2) );
  Eng -= 4361.42;
  Eng /= 939.57-sqrt( pow(1.0/_mom[1],2)+pow(105.66,2) )
    +fabs(1.0/_mom[1])*(1.0/sqrt( _dir[1][0]*_dir[1][0]+_dir[1][1]*_dir[1][1]+1 ));

  if ( Eng < 0 )//reconstruct as a free muon if -ve value calculated.
    Eng = sqrt( pow(1.0/_mom[1],2)+pow(105.66,2) );

  return Eng;
}

/********************************************/
double Golden_Analysis::quasi_nu_neut()
/********************************************/
{
  //Quasi formula: nu n scattering.
  double Eng;
  Eng = 939.57*sqrt( pow(1.0/_mom[1],2)+pow(105.66,2) );
  Eng -= 6802.61;
  Eng /= 938.27-sqrt( pow(1.0/_mom[1],2)+pow(105.66,2) )
    +fabs(1.0/_mom[1])*(1.0/sqrt( _dir[1][0]*_dir[1][0]+_dir[1][1]*_dir[1][1]+1));

  if ( Eng < 0 )//reconstruct as a free muon if -ve value calculated.
    Eng = sqrt( pow(1.0/_mom[1],2)+pow(105.66,2) );

  return Eng;
}

/********************************************/
double Golden_Analysis::calcQt(double sig)
/********************************************/
{
  double muVec[3], mag = 0, hadVec[3];
  //calculate muon unit direction vec.
  muVec[2] = 1.;
  muVec[0] = _dir[1][0];
  muVec[1] = _dir[1][1];
  
  for (int i=0;i<3;i++)
    mag += pow(muVec[i],2);

  mag = sqrt(mag);

  for (int j=0;j<3;j++)
    muVec[j] /= mag;
  
  //smear hadronic direction.
  double hadMag = sqrt( pow(_hadMom[0],2)+pow(_hadMom[1],2)+pow(_hadMom[2],2));
  double theta = acos( _hadMom[2]/hadMag );
  double cosThi = _hadMom[0]/(hadMag*sin( theta ));
  double sinThi = _hadMom[1]/(hadMag*sin( theta ));
  
  theta *= (180/3.1415926535897932);
  
  theta += sig;
  
  theta /= (180/3.1415926535897932);
  
  hadVec[0] = sin( theta )*cosThi;
  hadVec[1] = sin( theta )*sinThi;
  hadVec[2] = cos( theta );
  
  double ctheta = 0.;
  
  for (int k=0;k<3;k++)
    ctheta += muVec[k]*hadVec[k];
  
  return fabs(1.0/(_mom[1]*1000))*pow( sin( acos( ctheta ) ), 2 );
}

/********************************************/
bool Golden_Analysis::quad_and_transLong(int ID)
/********************************************/
{
  //Final cuts based on transverse vs longitudinal distance
  //and the refit of a parabola.
  int icand = 0;
  double lonLen, tranLen, qCharge;
  vector<double> XX, ZZ;

  while ( icand < _nhits && (int)XX.size() < _trHits[1] ){

    if ( _inCand[icand] ){
      XX.push_back( _xPos[icand] );
      ZZ.push_back( _zPos[icand] );
    }
    icand++;
  }
  
  //Longitudinal and transverse extent.
  lonLen = ZZ[ZZ.size()-1] - ZZ[0];
  tranLen = fabs( XX[XX.size()-1] - XX[0] );
  
  if ( tranLen/lonLen <= _translongCuts[0]
       -_translongCuts[1]*_trHits[1] ){
    _rcuts->Fill("Displacement",1); 
    if(ID==1) _rcuts_sig->Fill("Displacement",1); 
    else _rcuts_back->Fill("Displacement",1); 
    return false;
  }
  if ( lonLen <= _translongCuts[2] )
    if ( fabs(1.0/_mom[1]) > _translongCuts[3]*lonLen ){
      _rcuts->Fill("Displacement",1); 
      if(ID==1) _rcuts_sig->Fill("Displacement",1); 
      else _rcuts_back->Fill("Displacement",1); 
      return false;
    }
  
  _pcuts->Fill("Displacement",1); 
  if(ID==1) _pcuts_sig->Fill("Displacement",1); 
  else _pcuts_back->Fill("Displacement",1); 
 
  if(_quadFitCuts[0] != _quadFitCuts[1]){
    qCharge = fit_check( XX, ZZ );
    if(ID == 1) _qCharge->Fill(qCharge);
    else _qCharge_back->Fill(qCharge);
    if ( qCharge >= _quadFitCuts[0] && qCharge <= _quadFitCuts[1] ){
      _rcuts->Fill("Quadratic",1); 
      if(ID==1) _rcuts_sig->Fill("Quadratic",1); 
      else _rcuts_back->Fill("Quadratic",1); 
      return false;
    }
  }
  _pcuts->Fill("Quadratic",1); 
  if(ID==1) _pcuts_sig->Fill("Quadratic",1); 
  else _pcuts_back->Fill("Quadratic",1); 
  
  return true;
}


/********************************************/
double Golden_Analysis::fit_check(const vector<double>& XP, const vector<double>& ZP)
/********************************************/
{
  double Xfit[_trHits[1]], Zfit[_trHits[1]], Xerr[_trHits[1]], Zerr[_trHits[1]];
  
  for (unsigned int err = 0;err < XP.size();err++){
    
      Xfit[err] = XP[err];
      Zfit[err] = ZP[err];

      Xerr[err] = 10;
      Zerr[err] = 6;

  }
  
  TGraphErrors *gr1 = new TGraphErrors(_trHits[1], Zfit, Xfit, Zerr, Xerr);
  TF1 * func = new TF1("fit","pol2(0)",-3,3);
  gr1->Fit("fit", "QN");
  
  double c = func->GetParameter(2);
  double cE = func->GetParError(2);
  
  delete gr1;
  delete func;
  
  //+ve value corresponds to -ve charge
  if ( c/fabs(c) == -_charge[1] ) return fabs(cE/c);

  //return false;
  return -fabs(cE/c);//Ensures swaps have -ve and stays have +ve.

}

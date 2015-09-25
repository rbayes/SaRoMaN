#include <old_analysis.h>
//#include <TMath.h>

#include <TRandom3.h>

using namespace bhep;

//********************************************************************
old_analysis::old_analysis(bhep::prlevel level){
//********************************************************************

  m = bhep::messenger( level );

}

//********************************************************************
old_analysis::~old_analysis(){
//********************************************************************


}

//********************************************************************
bool old_analysis::initialize(const bhep::gstore& rstore,
			      const bhep::gstore& astore){
//********************************************************************

  m.message("+++ initializing stores +++",bhep::NORMAL);
  
  smear_store = rstore;
  
  cut_store = astore;
  
  //initialize random generator.
  long seed = (long)smear_store.fetch_dstore( "Gen_seed" );
  ranGen = TRandom3(seed);
  
  //initialize outputs.
  string filename = smear_store.fetch_sstore( "outFile" );
  outFile = new TFile( filename.c_str(), "recreate" );
  
  dataOut = new TTree("data", "Tree with mind smear analysis data");

  readParam();
  define_tree_branches();
  
  return true;
}

//********************************************************************
bool old_analysis::execute(const bhep::event& evt){
//********************************************************************
//perform smear and cuts.
  m.message("+++ Execute Function +++",bhep::NORMAL);
  
  bool ok;

  _parts.clear();
  get_event_properties( evt );

  get_particle_info();

  rec_Energy( evt );

  ok = cut_length();

  if ( ok )
    ok = kin_cuts();
  
  if ( ok ) _CC = 1;
  else _CC = 0;

  dataOut->Fill();

  return ok;
}

//********************************************************************
bool old_analysis::finalize(){
//********************************************************************

  //dataOut->Write();
  outFile->Write();
  outFile->Close();

  return true;
}

//********************************************************************
void old_analysis::readParam(){
//********************************************************************
//Get relevant parameters
  
  if ( smear_store.find_dstore("perC_sig") )
    _hadSmearSig1 = smear_store.fetch_dstore("perC_sig");

  if ( smear_store.find_dstore("const_sig") )
    _hadSmearSig2 = smear_store.fetch_dstore("const_sig");

  if ( smear_store.find_dstore("lamb_Fe") )
    _lambFe = smear_store.fetch_dstore("lamb_Fe");

  if ( smear_store.find_dstore("sqr_had_sig") )
    _had_ang_sig[0] = smear_store.fetch_dstore("sqr_had_sig");

  if ( smear_store.find_dstore("had_sig") )
    _had_ang_sig[1] = smear_store.fetch_dstore("had_sig");

  if ( smear_store.find_dstore("sample_dist") )
    _measRes = smear_store.fetch_dstore("sample_dist");

  if ( smear_store.find_dstore("x0") )
    _X0 = smear_store.fetch_dstore("x0");

  if ( smear_store.find_dstore("B_field") )
    _Bfield = smear_store.fetch_dstore("B_field");

  if ( cut_store.find_dstore("length_diff") )
    _minlenDiff = cut_store.fetch_dstore("length_diff");

  if ( cut_store.find_dstore("Qt") )
    _minQt = cut_store.fetch_dstore("Qt");

  if ( cut_store.find_dstore("mom_cut") )
    _momGrad = cut_store.fetch_dstore("mom_cut");
  
  // if ( cut_store.find_istore("sig_q") )
//     _sigCharge = cut_store.fetch_istore("sig_q");

}

//********************************************************************
void old_analysis::define_tree_branches(){
//********************************************************************
  
  dataOut->Branch("type", &_intType, "EventType/I");
  dataOut->Branch("ChargeCurrent", &_CC, "CC/B");
  dataOut->Branch("Eng", &_recEng, "recEng/D:truEng/D");
  dataOut->Branch("leadMom", &_leadp, "recp/D:trup/D");
  dataOut->Branch("Qt", &_Qt,"Qt/D");
  dataOut->Branch("deltaL", &_lendiff, "deltL/D");

}

//********************************************************************
void old_analysis::get_event_properties(const bhep::event& evt){
//********************************************************************
  
  _intType = get_int_type( evt.fetch_sproperty("IntType") );

  _parts = evt.true_particles();
  
  _recEng[1] = evt.fetch_dproperty("nuEnergy") / GeV;

}

//********************************************************************
int old_analysis::get_int_type(string intName){
//********************************************************************

  if ( intName == "CCQE" ) return 1; //CC quasi
  else if ( intName == "NCQE" ) return 2; //NC quasi
  else if ( intName == "1piRes" ) return 3; //single pion resonant
  else if ( intName == "miscRes" ) return 4; //other resonant
  else if ( intName == "CCDIS" ) return 5; //CC deep inelastic
  else if ( intName == "NCDIS" ) return 6; //NC deep inelastic
  else if ( intName == "CabQE" ) return 7; //Cabibo supressed quasi
  else if ( intName == "NCpi" ) return 8; //NC coherent/diff pi
  else if ( intName == "CCpi" ) return 9; //CC coherent/diff pi
  else if ( intName == "eEl" ) return 10; //elastic e- scattering
  else if ( intName == "muINVe" ) return 11; //inverse muon decay

}

//********************************************************************
void old_analysis::get_particle_info(){
//********************************************************************
//loop through particles to find the longest and second longest.
  vector<bhep::particle*>::iterator pIt;
  double Longest, no2Longest = 0;
  int nhits1, nhitsCurrent;
  double current;
  bool first = true;
  vector<const bhep::particle*> hijas;

  for (pIt = _parts.begin();pIt != _parts.end();pIt++){
    
    if ( (*pIt)->hits("tracking").size() != 0 )
      if ( first ){
	nhits1 = (int)(*pIt)->hits("tracking").size();
	Longest = (*pIt)->fetch_dproperty("length") / cm;
	_LpIt = pIt;
	hijas = (*pIt)->daughters();
	Longest += check_daughters( hijas );
	nhits1 += _dHits;
	first = false;
      } else {
	nhitsCurrent = (int)(*pIt)->hits("tracking").size();
	current = (*pIt)->fetch_dproperty("length") / cm;
	hijas = (*pIt)->daughters();
	current += check_daughters( hijas );
	nhitsCurrent += _dHits;
	
	//if (current > Longest){
	if (nhitsCurrent > nhits1){
	  no2Longest = Longest;
	  Longest = current;
	  _LpIt = pIt;
	} else if (current > no2Longest)
	  no2Longest = current;
      }
    hijas.clear();
  } 
  
  //Fill required properties. Still need to check units, need GeV.
  _leadL = Longest;
  _lendiff = Longest - no2Longest;
  _truMom = (*_LpIt)->p()/GeV;
  _leadp[1] = _truMom;
  _cosAng = ((*_LpIt)->px()/GeV)/_truMom;
  _nleadhits = (int)(*_LpIt)->hits("tracking").size();
  _leadLproj = (*_LpIt)->hits("tracking")[_nleadhits-1]
    - (*_LpIt)->hits("tracking")[0];
  _muEng = (*_LpIt)->e()/GeV;
  
}

//********************************************************************
double old_analysis::check_daughters(const vector<const bhep::particle*>& hijas){
//********************************************************************
//check for any decay daughters.
  double length = 0;
  _dHits = 0;
  vector<const bhep::particle*>::const_iterator dIt;

  for (dIt = hijas.begin();dIt != hijas.end();dIt++){

    if ( (*dIt)->fetch_sproperty("CreatorProcess")=="Decay" &&
	 (*dIt)->hits("tracking").size() != 0 ){
      if ( length == 0 ){
	length = (*dIt)->fetch_dproperty("length") / cm;
	_dHits = (int)(*dIt)->hits("tracking").size();
      } else
	if ( (*dIt)->fetch_dproperty("length") / cm > length ){
	  length = (*dIt)->fetch_dproperty("length") / cm;
	  _dHits = (int)(*dIt)->hits("tracking").size();
	}
    }

  }

  return length;
}

//********************************************************************
void old_analysis::rec_Energy(const bhep::event& evt){
//********************************************************************
//reconstruct neutrino energy.

  vector<bhep::particle*>::iterator pIt;

  bhep::Vector3D hadVec(0,0,0);

  bhep::Vector3D muVec = (*_LpIt)->p3() / GeV;
  
  _HadEng = 0;

  for (pIt = _parts.begin();pIt != _parts.end();pIt++)
    if ( (*pIt)->primary() )
      if ( pIt != _LpIt ){
	hadVec += (*pIt)->p3() / GeV;
	_HadEng += (*pIt)->e() / GeV;
      }

  double recHadE = rec_had_eng( _HadEng );

  double recMuE = rec_mu_eng();

  smear_ang( muVec, hadVec);

  _recEng[0] = recMuE + recHadE;

}

//********************************************************************
double old_analysis::rec_had_eng(double eng){
//********************************************************************
//Gaussian smear on total hadronic (had-lead in NC) energy
  double Esig = sqrt( pow(_hadSmearSig1*sqrt( eng ),2)
		      + pow (_hadSmearSig2*eng,2) );

  double EngRec = eng + ranGen.Gaus( 0, Esig);

  return EngRec;
}

//********************************************************************
double old_analysis::rec_mu_eng(){
//********************************************************************
//Gluckstern smear on lead particle momentum to reconstruct
//expected energy assuming muon
  
  double K, dKres, dKms;
  double res = _measRes;
  double beta = _truMom/_muEng;
  double sigMa, muP;
  
  K = (0.3*_Bfield)/(_cosAng * _truMom);

  dKres = ( res / pow( _leadLproj, 2) ) * sqrt( 720/( _nleadhits + 4 ));

  dKms = (0.016/(_leadL*_truMom*beta*pow(_cosAng,2)))*sqrt(_leadL/_X0);

  sigMa = sqrt( pow(dKres, 2) + pow(dKms, 2) );

  K += ranGen.Gaus( 0, sigMa);

  muP = (0.3*_Bfield)/(_cosAng*K);

  _leadp[0] = muP;

  return sqrt( pow(muP, 2) + pow(0.106, 2) );
}

//********************************************************************
void old_analysis::smear_ang(bhep::Vector3D& muVec,
			     bhep::Vector3D& hadVec){
//********************************************************************
//smear the polar angle.
  
  rec_mu_angle( muVec );

  rec_had_angle( hadVec );
  
  //_muhadang = muVec.angle( hadVec ); nolonger typedef from clhep so has less functionality
  _muhadang = acos( (muVec(0)*hadVec(0) + muVec(1)*hadVec(1) + muVec(2)*hadVec(2))
		    /(muVec.mag() * hadVec.mag()) );

}

//********************************************************************
void old_analysis::rec_mu_angle(bhep::Vector3D& muVec){
//********************************************************************

  double res = _measRes;
  double beta = _truMom/_muEng;
  double cc = 300000000;
  //get polar angles.
  double theta = acos( muVec[2]/_truMom );
  double cosThi = muVec[0]/(_truMom*sin( theta ));
  double sinThi = muVec[1]/(_truMom*sin( theta ));
  double dTHres, dTHms, dTH, Nbar = 0;
  
  for (float jj = _nleadhits;jj > 0;jj--)
    Nbar += pow(jj, 2);

  Nbar = sqrt( Nbar );

  theta *= (180 / Pi);//to degrees

  dTHres = res/( _lambFe * Nbar );
  
  dTHms = (0.0136/(beta*cc*_truMom)) * sqrt( (_lambFe*_nleadhits)/_X0 )
    * (1 + 0.038*log( (_lambFe*_nleadhits)/_X0 ));
  
  dTH = sqrt( pow(dTHres, 2) + pow(dTHms, 2) );

  theta += ranGen.Gaus( 0, dTH);
  theta /= (180 / Pi);//back to rad.
  
  muVec[0] = _leadp[0] * sin( theta ) * cosThi;
  muVec[1] = _leadp[0] * sin( theta ) * sinThi;
  muVec[2] = _leadp[0] * cos( theta );

}

//********************************************************************
void old_analysis::rec_had_angle(bhep::Vector3D& hadVec){
//********************************************************************

  double hadMomTru = hadVec.mag();
  double theta = acos( hadVec[2]/hadMomTru );
  double cosThi = hadVec[0]/(hadMomTru*sin( theta ));
  double sinThi = hadVec[1]/(hadMomTru*sin( theta ));
  double dTH;

  theta *= (180 / Pi);//to degrees

  dTH = sqrt( pow(_had_ang_sig[0]/sqrt( _HadEng ), 2)
	      + pow(_had_ang_sig[1]/_HadEng, 2) );

  theta += ranGen.Gaus( 0, dTH);
  theta /= (180 / Pi);//back to rad.
  
  hadVec[0] = hadMomTru * sin( theta ) * cosThi;
  hadVec[1] = hadMomTru * sin( theta ) * sinThi;
  hadVec[2] = hadMomTru * cos( theta );

}

//********************************************************************
bool old_analysis::cut_length(){
//********************************************************************

  if ( _lendiff < _minlenDiff ){

    _CC = 0;
    _Qt = -1;//not calculated so set to negative value in tree

    return false;

  }

  return true;
}

//********************************************************************
bool old_analysis::kin_cuts(){
//********************************************************************

  bool ok;
  
  _Qt = _leadp[0] * pow( sin( _muhadang ), 2);
  double momCut = _momGrad * _recEng[0];
  cout << "cut levels: "<<_Qt << ", "<<_minQt<<", mom: "<<_leadp[0]<<", "<<momCut<< endl;
  if ( _recEng[0] < 7 ) ok = true; //assumming GeV.
  else if ( _Qt < _minQt || _leadp[0] < momCut ) ok = false;
  else ok = true;

  return ok;
}

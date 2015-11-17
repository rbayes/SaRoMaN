
#include <root2dst.h>
#include <digi/WriteUtil.h>

#include <bhep/bhep_svc.h>
//#include <CLHEP/Random/RandGauss.h>

using namespace std;
using namespace bhep;


//*************************************************************
root2dst::root2dst(bhep::prlevel vlevel, bhep::gstore* store){
//*************************************************************
  
  level = vlevel;
  
  m = bhep::messenger(level);

  dataIn = NULL;

  if ( store != NULL ){
    _construct = new gdml_hit_constructor( *store );
    _doVox = true;
  } else _doVox = false;
  
  m.message("+++ root2dst Constructor +++",bhep::VERBOSE); 
  
}


//*************************************************************
bool root2dst::initialize(double *res, long seed,
			  TTree *InPutTree, TString OutFileName) {
//*************************************************************
    
  m.message("+++ root2dst init  function ++++",bhep::NORMAL);
  
  if (InPutTree != NULL){
    //outgz.open((string)OutFileName);
    WriteUtil::OpenOutputDst( (string)OutFileName );
    
    dataIn = InPutTree;
    dataIn->SetBranchStatus("*",0);

  }

  //long seed = 178263094;
  //ranGen = RanluxEngine(seed, 4);
  //cout << res[0] << "," << res[1] << endl;
  sigMa = res[0];
  sigMaE = res[1];

  nevt=0;
  
  return true;
}


//*************************************************************
bool root2dst::execute(){
//*************************************************************
    
  /*
    Take a ROOT tree event and convert it into bhep DST event
   */

  m.message("+++ root2dst execute function ++++",bhep::VERBOSE);
  cout << nevt<<endl;
  createEvent();
  
  bool ok = make_particles();

  if ( !ok ) std::cout << "Bad Particles in event" << std::endl;

  //outgz.write(nuEvent, nevt);
  bhep::bhep_svc::instance()->get_writer_root().write( nuEvent, nevt );
  nevt++;

  nuEvent.clear();
    
  return true;

}



//*************************************************************
bool root2dst::finalize() {
//*************************************************************
  
  if (dataIn != NULL){
    //outgz.close();
    WriteUtil::CloseOutputDst();
    delete dataIn;
  }

  m.message("+++ root2dst finalize function ++++",bhep::NORMAL);
    
  m.message("++ Number of analyzed events: ",nevt,bhep::NORMAL);

  return true;
}


//*************************************************************
void root2dst::createEvent() {
//*************************************************************

/* Defines an event and adds the event specific information. */

  Int_t integs[2][2];
  Float_t floas[6][3];

  nuEvent.set_event_number((int)nevt);

  dataIn->SetBranchStatus("Random",1); dataIn->SetBranchAddress("Random",&integs[0]);
  dataIn->SetBranchStatus("Nutype",1); dataIn->SetBranchAddress("Nutype",&integs[1]);
  dataIn->SetBranchStatus("Enu",1); dataIn->SetBranchAddress("Enu",&floas[0]);
  dataIn->SetBranchStatus("Xf",1); dataIn->SetBranchAddress("Xf",&floas[1]);
  dataIn->SetBranchStatus("Yf",1); dataIn->SetBranchAddress("Yf",&floas[2]);
  dataIn->SetBranchStatus("Q2f",1); dataIn->SetBranchAddress("Q2f",&floas[3]);
  dataIn->SetBranchStatus("Wf",1); dataIn->SetBranchAddress("Wf",&floas[4]);
  dataIn->SetBranchStatus("Vertex",1); dataIn->SetBranchAddress("Vertex",&floas[5]);
  
  dataIn->GetEntry((Int_t)nevt);

  nuEvent.set_vertex(floas[5][0] * cm, floas[5][1] * cm, floas[5][2] * cm);
  
  //Add aditional properties
  nuEvent.add_property("RandomSeed1", integs[0][0] );
  nuEvent.add_property("RandomSeed2", integs[0][1] );
  nuEvent.add_property("Nutype", integs[1][0] );
  nuEvent.add_property("nuEnergy", (double)floas[0][0] * GeV );
  nuEvent.add_property("Xf", (double)floas[1][0] );
  nuEvent.add_property("Yf", (double)floas[2][0] );
  nuEvent.add_property("Q2f", (double)floas[3][0] );
  nuEvent.add_property("Wf", (double)floas[4][0] );
  
  dataIn->SetBranchStatus("*", 0);

  cout <<"Event Defined" << endl;
}

//***************************************************************
bool root2dst::make_particles() {
//***************************************************************

/* Function to add relevant particles to event */
  particle *lead, *hadron, *digPar;
  string detect = "MIND";

  lead = define_lead_particle();
  
  hadron = define_hadron();

  vector<hit*> muHit, hadHit;
  bool ok = hits_fromFile(muHit, hadHit);
  vector<particle*> all_parts;
  
  if (ok) {

    if (muHit.size() != 0)
      for (Int_t iHit=0;iHit<(Int_t)muHit.size();iHit++){
	if ( lead->name() == "mu-" || lead->name() == "mu+" ){
	  muHit[iHit]->set_mother_particle(*lead);
	  lead->add_hit(detect, muHit[iHit]);
	} else {
	  muHit[iHit]->set_mother_particle(*hadron);
	  hadron->add_hit(detect, muHit[iHit]);
	}
      }
    if (hadHit.size() != 0)
      for (Int_t jHit=0;jHit<(Int_t)hadHit.size();jHit++){
	hadHit[jHit]->set_mother_particle(*hadron);
	hadron->add_hit(detect, hadHit[jHit]);
      }
  
    all_parts.push_back( lead );
    all_parts.push_back( hadron );

    digPar = create_digital_representation(all_parts);
    
    nuEvent.add_true_particle(lead);
    nuEvent.add_true_particle(hadron);
    nuEvent.add_digi_particle(digPar);

  }
  else cout << "No. hits in event" << endl;

  cout << "All particles defined and appended to event" << endl;
  return ok;
}

//***************************************************************
particle* root2dst::define_lead_particle() {
//***************************************************************

/* Function to define the properties of the lead particle from
   the event. e.g. the mu- from a numucc interaction */
  
  Int_t ID;
  Float_t moment[5];
  
  ptype pT = TRUTH;
  string pName;

  dataIn->SetBranchStatus("Idlead", 1); dataIn->SetBranchAddress("Idlead", &ID);
  dataIn->SetBranchStatus("Plead", 1); dataIn->SetBranchAddress("Plead", &moment);
  dataIn->SetBranchStatus("Nhits", 1); dataIn->SetBranchAddress("Nhits", &NmuHits);

  dataIn->GetEntry((Int_t)nevt);
  
  Point3D pos = nuEvent.vertex();
  //multiply by 1000 because ntuple in GeV and DST expects MeV.
  Vector3D mom(moment[0] * GeV, moment[1] * GeV, moment[2] * GeV);
  
  ray *leadstart = new ray(pos, mom);
  
  switch (ID) {

  case 13:
    pName = "mu-";
    break;
  case -13:
    pName = "mu+";
    break;
  case 14:
    pName = "nu_mu";
    break;
  case -14:
    pName = "anti_nu_mu";
    break;
  case 11:
    pName = "e-";
    break;
  case -11:
    pName = "e+";
    break;
    
  }
  
  particle *lead = new particle(pT, pName, *leadstart);

  dataIn->SetBranchStatus("*", 0);

  delete leadstart;
  
  return lead;
  
}

//**************************************************************
particle* root2dst::define_hadron() {
//**************************************************************

/* Add particle to represent the hadronic shower associated 
   with the event */
  
  ptype pT = TRUTH;
  string pName = "Hadronic_vector";

  Int_t nHad;
  Float_t moment[5];
  Float_t hadEng;

  dataIn->SetBranchStatus("Nhad", 1); dataIn->SetBranchAddress("Nhad", &nHad);
  dataIn->SetBranchStatus("Phadg", 1); dataIn->SetBranchAddress("Phadg", &moment);
  dataIn->SetBranchStatus("Ehadg", 1); dataIn->SetBranchAddress("Ehadg", &hadEng);
  dataIn->SetBranchStatus("Nhhits", 1); dataIn->SetBranchAddress("Nhhits", &NhadHits);

  dataIn->GetEntry((Int_t)nevt);
  
  Point3D pos = nuEvent.vertex();
  
  Vector3D mom(moment[0] * GeV, moment[1] * GeV, moment[2] * GeV);
  
  ray *hadstart = new ray(pos, mom);
  
  particle *had = new particle(pT, pName, *hadstart);
  cout << "Particle Defined" << endl;
  had->add_property("No_Hadrons_in_Shower", nHad );
  had->add_property("HadE", (double)hadEng );

  dataIn->SetBranchStatus("*", 0);

  delete hadstart;

  //cout << "Hadronic shower  vector defined" << endl;

  return had;
}


//**************************************************************************
bool root2dst::hits_fromFile(vector<hit*>& muHit, vector<hit*>& hadHit) {
//**************************************************************************

/* Add vector of hit locations associated with particle par */

  const Int_t maxHits = 400;
  string detect = "MIND";

  //Int_t NmuHits, NhadHits;
  Float_t muHitVec[3][maxHits], hadHitVec[3][maxHits];
  Float_t muHitE[maxHits], hadHitE[maxHits];
  
  if (NmuHits==0 && NhadHits==0) return false;

  //dataIn->SetBranchStatus("Nhits", 1); dataIn->SetBranchAddress("Nhits", &NmuHits);
  dataIn->SetBranchStatus("Xhit", 1); dataIn->SetBranchAddress("Xhit", &muHitVec[0]);
  dataIn->SetBranchStatus("Yhit", 1); dataIn->SetBranchAddress("Yhit", &muHitVec[1]);
  dataIn->SetBranchStatus("Zhit", 1); dataIn->SetBranchAddress("Zhit", &muHitVec[2]);
  dataIn->SetBranchStatus("Ehit", 1); dataIn->SetBranchAddress("Ehit", &muHitE);
  //dataIn->SetBranchStatus("Nhhits", 1); dataIn->SetBranchAddress("Nhhits", &NhadHits);
  dataIn->SetBranchStatus("Xhhit", 1); dataIn->SetBranchAddress("Xhhit", &hadHitVec[0]);
  dataIn->SetBranchStatus("Yhhit", 1); dataIn->SetBranchAddress("Yhhit", &hadHitVec[1]);
  dataIn->SetBranchStatus("Zhhit", 1); dataIn->SetBranchAddress("Zhhit", &hadHitVec[2]);
  dataIn->SetBranchStatus("Ehhit", 1); dataIn->SetBranchAddress("Ehhit", &hadHitE);

  dataIn->GetEntry((Int_t)nevt);

  if (NmuHits!=0)
    for (Int_t iHit = 0;iHit < NmuHits;iHit++){
      
      Point3D hitPos(muHitVec[0][iHit] * cm, muHitVec[1][iHit] * cm, muHitVec[2][iHit] * cm);
      
      hit *xyz = new hit(detect);
      xyz->set_point(hitPos);
      //xyz->add_data("E_dep", bhep::to_string( muHitE[iHit] ));
      xyz->add_property("EnergyDep", (double)muHitE[iHit] );

      muHit.push_back(xyz);
      
    }

  if (NhadHits!=0)
    for (Int_t jHit = 0;jHit < NhadHits;jHit++){
      
      Point3D hitPos(hadHitVec[0][jHit] * cm, hadHitVec[1][jHit] * cm, hadHitVec[2][jHit] * cm);
      
      hit *xyz = new hit(detect);
      xyz->set_point(hitPos);
      //xyz->add_data("E_dep", bhep::to_string( hadHitE[jHit] ));
      xyz->add_property("EnergyDep", (double)hadHitE[jHit] );

      hadHit.push_back(xyz);
      
    }

  dataIn->SetBranchStatus("*", 0);

  return true;
}


//***************************************************************
particle* root2dst::create_digital_representation(const vector<particle*>& tru_parts) {
//***************************************************************

/*Makes a digital particle with only the hits, to be used in
  the fitting algorithm */


  ptype pT = DIGI;
  string detect = "tracking";
  string edep = "EnergyDep";//"E_dep";

  particle *hitMap = new particle(pT, "unknown");

  vector<bhep::hit*> part_hits, temp_hit;
  int moth_prim;
  double X, Y, Z;
  double truE, recE, E_sig;

  if ( _doVox ){

    vector<bhep::hit*> vox;
    vector<bhep::hit*>::iterator voxIt;
    
    for (int iPart = 0;iPart < (int)tru_parts.size();iPart++){
      temp_hit = tru_parts[iPart]->hits(detect);
      
      if ( part_hits.size() == 0 && temp_hit.size() != 0 )
	part_hits = temp_hit;
      else if ( temp_hit.size() != 0 )
	part_hits.insert( part_hits.end(), temp_hit.begin(), temp_hit.end() );

      temp_hit.clear();
    }
    
    _construct->execute( part_hits, vox );
    
    for (voxIt = vox.begin();voxIt != vox.end();voxIt++)
      hitMap->add_hit( detect , (*voxIt) );

  } else {
    for (Int_t iPart = 0;iPart < (Int_t)tru_parts.size();iPart++) {
      
      part_hits = tru_parts[iPart]->hits( detect );
      
      if (part_hits.size()!=0)
	for (Int_t iHit = 0;iHit < (Int_t)part_hits.size();iHit++) {
	  
	  //truE = bhep::double_from_string( part_hits[iHit]->data("E_dep") );
	  truE = part_hits[iHit]->ddata( edep );
	  
	  //E_sig = sigMaE * truE;
	  //if ( truE > 0.0009 ){
	  
	  recE = truE;// + RandGauss::shoot(&ranGen, 0, E_sig);
	  
	  X = part_hits[iHit]->x().x()/cm;// + RandGauss::shoot(&ranGen, 0, sigMa);
	  Y = part_hits[iHit]->x().y()/cm;// + RandGauss::shoot(&ranGen, 0, sigMa);
	  Z = part_hits[iHit]->x().z()/cm;
	  
	  Point3D hitPos(X * cm,Y * cm,Z * cm);
	  hit* digHit = new hit( detect );
	  digHit->set_point(hitPos);
	  //digHit->add_data("E_dep", bhep::to_string( recE ));
	  digHit->add_property( edep , recE);
	  digHit->add_property( "true_moth", tru_parts[iPart]->name() );
	  if ( part_hits[iHit]->mother_particle()
	       .fetch_sproperty("CreatorProcess")=="none" )
	    moth_prim = 1;
	  else
	    moth_prim = 0;
	  digHit->add_property( "moth_prim", moth_prim );
	  //cout << tru_parts[iPart]->name() << endl;
	  //digHit->set_mother_particle( *tru_parts[iPart] );
	  hitMap->add_hit( detect , digHit);
	  //}
	}
      part_hits.clear();
      
    }
  }
  
  return hitMap;
}

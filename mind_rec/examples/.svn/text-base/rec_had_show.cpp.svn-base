/***********************************************************
 *                                                         *
 * Example to study reconstruction of hadron shower from   *
 * a neutrino DIS event.                                   *
 *                                                         *
 * June 2009.                                              *
 *                                                         *
 ***********************************************************/

#include <mind/shower_rec.h>
#include <mind/shower.h>

//#include <bhep/EventManager2.h>
#include <bhep/gstore.h>
#include <bhep/event.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>

#include <recpack/Measurement.h>
#include <recpack/stc_tools.h>

#include <TFile.h>
#include <TTree.h>

void get_event_data(bhep::event& evt, measurement_vector& meas, EVector& vert, double* mom);
Measurement* getMeasurement(bhep::hit& hit);

int main(int argc, char *argv[]){
  
  if ( argc < 3 ){

    cout << "execute as ./rec_had_show <param file> <no. events>" <<endl;
    exit(1);

  }
  
  bool succ;
  double truP[3], recDir[3];
  
  string param_file = argv[1];

  int nEvent = atoi( argv[2] );
  
  bhep::gstore data_store;

  bhep::sreader reader(data_store);
  reader.file(param_file);
  reader.group("DATA"); reader.read();
  bhep::sreader reader2(data_store);
  reader2.file(param_file);
  reader2.group("ANA");reader2.read();
  
  //EventManager2 *eman = new EventManager2(data_store,bhep::MUTE);

  //eman->initialize();
  bhep::reader_root inDst;
  inDst.open( data_store.fetch_sstore("idst_file") );
  
  shower_rec* recon = new shower_rec();
  
  recon->initialize( data_store.fetch_dstore("tol"), data_store.fetch_dstore("res") );
  
  TFile *output = new TFile("test_out.root", "recreate");
  TTree *data = new TTree("t1","had. rec. info.");
  
  data->Branch("Success",&succ,"succ/B");
  data->Branch("TrueMom",&truP,"truP[3]/D");
  data->Branch("RecDirection",&recDir,"dir[3]/D");
  
  measurement_vector meas;
  EVector vert(3,0);
  
  for (int i = 0;i < nEvent;i++){
    
    //succ = eman->status();

    //if ( !succ ) return -1;
    //bhep::event& e = eman->read();
    bhep::event& e = inDst.read_event( i );
    std::cout << "...got event " << i <<std::endl;
    get_event_data( e, meas, vert, truP );
    
    if ( meas.size() != 0 ){

      shower *had_show = new shower( meas, vert );
      
      succ = recon->execute( *had_show );
      
      recDir[0] = had_show->get_direction()[0];
      recDir[1] = had_show->get_direction()[1];
      recDir[2] = had_show->get_direction()[2];
      
      data->Fill();
      
      delete had_show;
      stc_tools::destroy( meas );
    }
    
  }
  
  output->Write();
  output->Close();
  
  recon->finalize();
  
  //eman->finalize();
  inDst.close();
  
  return 0;
}

void get_event_data(bhep::event& evt, measurement_vector& meas, EVector& vert, double* mom)
{
  vert[0] = evt.vertex()[0];
  vert[1] = evt.vertex()[1];
  vert[2] = evt.vertex()[2];

  //make the measurement vector.
  vector<bhep::particle*> parts = evt.true_particles();
  vector<bhep::particle*>::iterator pIt;

  //in this version use the truth. the combined hadron particle from G3
  for (pIt = parts.begin();pIt != parts.end();pIt++)
    {
      if ( (*pIt)->name().compare("Hadronic_vector") == 0 ){

	const vector<bhep::hit*>& hits = (*pIt)->hits("MIND");
	mom[0] = (*pIt)->px();
	mom[1] = (*pIt)->py();
	mom[2] = (*pIt)->pz();

	for (size_t j=0; j< hits.size(); j++){

	  Measurement* mnt = getMeasurement( *hits[j] );

	  meas.push_back( mnt );

	}
      }
    }

}

Measurement* getMeasurement(bhep::hit& hit)
{
  EVector hit_pos(2,0);
  EVector meas_pos(3,0);
  
  meas_pos[0] = hit.x()[0];
  meas_pos[1] = hit.x()[1];
  meas_pos[2] = hit.x()[2];
  
  hit_pos[0] = meas_pos[0];
  hit_pos[1] = meas_pos[1];

  //covariance and meastype hardwired for now.
  EMatrix cov(2,2,0);
  cov[0][0] = 1.; cov[1][1] = 1.;
  string meastype = "xy";

  Measurement* me = new Measurement();
  me->set_name(meastype);
  me->set_hv(HyperVector(hit_pos,cov));
  me->set_name("volume", "Detector");
  me->set_position( meas_pos );

  //Add the hit energy deposit as a key to the Measurement.
  const dict::Key Edep = "E_dep";
  const dict::Key EdepVal = hit.data("E_dep");
  me->set_name(Edep, EdepVal);

  return me;
}

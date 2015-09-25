#include <mind/fitter.h>
#include <mind/MINDplotter.h>

#include <bhep/gstore.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>

#include <recpack/Ring.h>
#include <recpack/stc_tools.h>

using namespace std;

/*


 */

int main(int argc, char* argv[]){
    
  //---- set parameters file ----//
    
  string param_file="param/fit_cats.param";
  
  int nevents=1;
  bool fitOk;
  if (argc==2) param_file = argv[1];
  
  else if(argc==3) {param_file = argv[1];nevents = atoi(argv[2]);}
  
  else{cout<<
	 "\n++ Usage: ./fit_tracks param_file number_of_events\n"<<endl;
  exit(1);}
  
  // generate stores for analysis parameters 
  bhep::gstore data_store,ana_store;//,run_store;
  
  // read analysis parameters 
  bhep::sreader reader(ana_store);
  reader.file(param_file);
  reader.group("ANA");reader.read();
  bhep::sreader reader2(ana_store);
  reader2.file(param_file); 
  reader2.group("RUN");reader2.read();
  //
  
  // // read run properties
//   bhep::sreader preader(run_store);
//   preader.file(param_file);
//   preader.group("RUN");
//   preader.read();
//   //
  
  // read input/ouput data
  bhep::sreader data_reader(data_store);
  data_reader.file(param_file);
  data_reader.group("DATA");
  data_reader.read();
  //
  
  bhep::reader_root inDst;
  
  // fitter* fit = new fitter(ana_store,bhep::MUTE);
  fitter fit(ana_store,bhep::MUTE);
  
  //  MINDplotter* plot = new MINDplotter();
  MINDplotter plot = MINDplotter();
  
  //catchOk = fit->initialize();
  fit.Initialize();
  
  bool patR = ana_store.fetch_istore("patRec");
  bool clusts = ana_store.fetch_istore("do_clust");
  
  //catchOk = plot->initialize(run_store.fetch_sstore("out_file"),bhep::MUTE);
  plot.Initialize(ana_store.fetch_sstore("out_file"), patR, clusts, bhep::MUTE);

  /// changed to get CA output
  //  plot.initialize(ana_store.fetch_sstore("out_file"), ana_store.fetch_sstore("CA_file"), patR, clusts, bhep::MUTE);
  
  vector<string> input_data = data_store.fetch_svstore("idst_files");
  
  //Counters for event loops;
  int i;
  int evt_read = 0;
  string filetest;
  //
  
  for (unsigned int ifile = 0;ifile < input_data.size();ifile++){
    
    inDst.open( input_data[ifile] );
    i = 0;
    
    //for(int i=0; i < nevents; i++) {
    while ( !inDst.eof(i) && evt_read < nevents ) {
      
      // if (i%100==0) cout<< "Number of events read "<<evt_read<<endl;
      cout << "Event: " << i << endl;
      bhep::event& e = inDst.read_event( i );
      if ( e.find_sproperty("IntType") ){//Protects against corrupt events 
	//caused by G4_out being spread over more than one file.
	// loop over particles
	vector<bhep::particle*> parts = e.digi_particles();
	
	//Relevant only when building likelihood tree.
	fit.set_int_type( e.fetch_sproperty("IntType") );
	//
	// cout <<"There are " << parts.size() << " digis in event "
	// 	     << e.event_number() <<endl;
	
	if (parts.size() != 0) {
	  for (size_t part=0; part<parts.size();part++){
	    
	    if (parts[part]->name()=="void") continue;
	    
	    fitOk = fit.Execute(*parts[part],e.event_number());
	    
	    ///for single track plotter
	    // plot.execute(fit, e, fitOk);
	    
	    ///for multiple track plotter
	    plot.Execute(fit, e);
	    
	  }
	}
	
	parts.clear();
	e.clear();
	
	//i++;
	evt_read++;
      }
      i++;
    }
    
    inDst.close();
    
  }
  
  fit.Finalize();
  
  plot.Finalize();
  
  
  std::cout<<"Reconstruction Complete\n";
  
  return 0;
  
}

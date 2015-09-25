/************************************************************
 *                                                          *
 * Program to control the smear and cut based analysis      *
 * on mind Geant4 bhep output.                              *
 *                                                          *
 * Execute as:                                              *
 *        ./smear_analysis <param_file_name> <noEvents>     *
 *                                                          *
 * Authors: A. Laing, A. Cervera; February 2009.            *
 *                                                          *
 ************************************************************/

#include <digi/old_analysis.h>
#include <digi/WriteUtil.h>

//#include <bhep/EventManager2.h>
#include <bhep/gstore.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
  
  string param_file;
  int nEvents;

  if (argc == 3) {
    
    param_file = argv[1];
    nEvents = atoi( argv[2] );

  } else {

    std::cout << "Incorrect execution" << std::endl;

    return -1;

  }
  
  bhep::gstore run_store, ana_store, data_store;

  bhep::sreader reader1( run_store );
  reader1.file( param_file );
  reader1.group( "RUN" );
  reader1.read();

  bhep::sreader reader2( ana_store );
  reader2.file( param_file );
  reader2.group( "ANA" );
  reader2.read();

  bhep::sreader reader3( data_store );
  reader3.file( param_file );
  reader3.group( "DATA" );
  reader3.read();

  //EventManager2* eman = new EventManager2( data_store, bhep::VERBOSE);
  bhep::reader_root inDst;
  inDst.open( data_store.fetch_sstore("idst_file") );
  
  old_analysis* ays = new old_analysis( bhep::NORMAL );
  
  //bool shite = eman->initialize();
  
  ays->initialize( run_store, ana_store);
  
  // if ( nEvents > (int)eman->getNumEventsFile() )
  //     nEvents = (int)eman->getNumEventsFile();
  
  for (int iEv = 0;iEv < nEvents;iEv++) {
    
    //bhep::event& e = eman->read();
    bhep::event& e = inDst.read_event( iEv );
    
    //if (iEv==0) continue;
    
    ays->execute( e );
    
  }
  
  ays->finalize();
  //eman->finalize();
  inDst.close();
  
  return 0;
}

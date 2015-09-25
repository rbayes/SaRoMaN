#include "Golden_Analysis.h"

#include <bhep/gstore.h>
#include <bhep/sreader.h>

int main(int argc, char* argv[]) {

  std::string param_file;
  std::string outFileName;

  if ( argc == 4 ) { param_file = argv[1]; outFileName = argv[2]; }
  else {

    std::cout << "Execute as: ./run_analysis <parameter file> <outfile> <seed>" << std::endl;

    return -1;
  }

  bhep::gstore run_store;

  //Read run and cut parameters.
  bhep::sreader reader1(run_store);
  reader1.file(param_file);
  reader1.group("RUN");
  reader1.read();
  bhep::sreader reader2(run_store);
  reader2.file(param_file);
  reader2.group("CUT");
  reader2.read();
  //Add the seed to the store.
  run_store.store("outFile", outFileName );
  run_store.store("seed", strtod( argv[3], NULL ) );
  //
  
  Golden_Analysis run1( run_store );
  
  run1.execute();
  
  run1.finalize();
  
  return 0;
}

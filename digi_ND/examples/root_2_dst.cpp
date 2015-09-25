
#include <digi/root2dst.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;


/*****************************************************************
 * Root tree to bHEP dst converter. Takes in a root file reads   *
 * the information in the trees and creates a gz file with the   *
 * relevant events as objects.                                   *
 *                                                               *
 * Execution via calls of the form:                              *
 *  ./root_2_dst <Input Filename> <Gaussian smear sig/cm>        *
 *  with option of specifying <No. events>                       *
 *                                                               *
 * Authors: Andrew Laing, Pau Novella.                           *
 *****************************************************************/

int main(int argc, char* argv[]){
    
  if (argc<5){
    
    cout << "Execution requires 2 or 3 arguments." << endl;
    cout << "Call with ./root_2_dst <InFile> <Gaus sigma, in cm> <E res> <rndm seed>" << endl;
    cout << "and optional  <No. events of interest>" << endl;

    return -1;
  }
  
  Char_t *inFileName;
  double smearRes[2];
  long rndmSeed;
  Int_t nEvents;

  inFileName = argv[1];

  TString outFileName = inFileName;
  outFileName.Replace(outFileName.Length()-4, 4, "dst.root");

  //Retreive Tree from file.
  TTree *Data = NULL;

  TFile In(inFileName);
  In.GetObject("h10", Data);

  //Get Gaussian resulution for smear.
  smearRes[0] = atof(argv[2]);

  //Get Energy smear.
  smearRes[1] = atof(argv[3]);

  //Get seed value for random engine
  rndmSeed = (long)atof(argv[4]);

  //Set number of events to be read.
  if (argc==6) nEvents = atoi(argv[5]);
  else nEvents = (Int_t)Data->GetEntries();

  if (nEvents>(Int_t)Data->GetEntries()){
    cout << "Input no. events greater than available in file.\n"
	 << "Replacing with max entries" << endl;
    nEvents = (Int_t)Data->GetEntries();
  }
  
  bhep::prlevel c = bhep::VERBOSE;
    
  root2dst* cvt = new root2dst(c);
    
  cvt->initialize(smearRes, rndmSeed, Data, outFileName);
  
  cout << "Starting Loop over events" << endl;
  for(int i=1; i < nEvents+1; i++) {
    
    if (i%100==0) cout<< "Number of events read "<<i<<endl;
    
    cvt->execute();
    
  }
  
  cvt->finalize();
  
  return 0;
  
}

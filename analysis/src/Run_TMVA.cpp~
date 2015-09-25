#include <sstream>
#include "TRandom3.h"

#include "Run_TMVA.h"

Run_TMVA::Run_TMVA(const bhep::gstore& store)
/*: histLk(0), histLkD(0), histLkPCA(0), histLkKDE(0), histLkMIX(0), histPD(0), histPDD(0),
    histPDPCA(0), histPDEFoam(0), histPDEFoamErr(0), histPDEFoamSig(0), histKNN(0), histHm(0),
    histFi(0), histFiG(0), histFiB(0), histLD(0), histNn(0),histNnbfgs(0),histNnbnn(0),
    histNnC(0), histNnT(0), histBdt(0), histBdtG(0), histBdtD(0), histRf(0), histSVMG(0),
    histSVMP(0), histSVML(0), histFDAMT(0), histFDAGA(0), histCat(0), histPBdt(0),
    recEffKNN(0), recEffLD(0), recEffBdt(0), recBackKNN(0), recBackLD(0), recBackBdt(0),
    recEffKNNe(0), recEffLDe(0), recEffBdte(0), recBackKNNe(0), recBackLDe(0), recBackBdte(0),
    alltracks2(0)
*/  
{
  
  // Initialize Classification step of analysis.
  
  // The output file name
  std::string outName = store.fetch_sstore("outFile");
  _OutFile = new TFile( outName.c_str(), "RECREATE" );
  // The input file bases for the training of the trees
  _filebaseTree  = store.fetch_sstore("TreeFileBase");
  _firstT     = store.fetch_istore("firstTree");
  _lastT      = store.fetch_istore("lastTree");
  int seed = store.fetch_istore("seed");


  rand = new TRandom3(seed);
  set_MVA_methods(store);
  _reader = new TMVA::Reader( "!Color:!Silent" );
  // Set the fixed cuts based on the config file
  define_cuts(store);
  define_variables(store);
  define_histograms(store);

}


void Run_TMVA::set_MVA_methods(const bhep::gstore& store)
{
  // By default only the k-nearest neighbour method will be used
  //std::map<std::string,int> Use;

  // --- Cut optimisation
  Use["Cuts"]            = 0;
  Use["CutsD"]           = 0;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 0;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 0;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 0; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 0; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 0;
  // 
  // --- Boosted Decision Trees
  Use["BDT"]             = 0; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 0;
  
  // If another method is to be used it can be read from the configuration file 
  if( store.find_svstore("MethodList")){
    _myMethodList = store.fetch_svstore("MethodList");
    
    for (UInt_t i=0; i<_myMethodList.size(); i++) {
      std::string regMethod(_myMethodList[i]);
      
      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
	return;
      }
      Use[regMethod] = 1;
      std::cout<<"Will use the "<<regMethod<<" method in analysis."<<std::endl;
    }
  }
}

void Run_TMVA::define_cuts(const bhep::gstore& store){
  
  _Vtype      = store.find_istore("vtype") ?
    store.fetch_istore("vtype"): 0;
  _doVary     = store.find_istore("doVary") ?
    store.fetch_istore("doVary"): 0;
  _beamCharge = store.fetch_istore("bChar");
  
  _maxE       = store.fetch_dstore("truEmax");
  _maxOver    = store.fetch_dstore("maxOver");

  _minNodes   = store.fetch_dstore("nodeprop");

  if(store.find_vstore("MVAcut"))
    _MVAcut     = store.fetch_vstore("MVAcut");

  _smearRoot = store.find_dstore("rootP") ?
    store.fetch_dstore("rootP") : 0.55;
  _smearConst = store.find_dstore("constP") ?
    store.fetch_dstore("constP") : 0.03;

  _smearPr = store.find_dstore("PrBlur") ?
    store.fetch_dstore("PrBlur") : 0;
  _smearPp = store.find_dstore("PpBlur") ?
    store.fetch_dstore("PpBlur") : 0;
  _smearPk = store.find_dstore("PkBlur") ?
    store.fetch_dstore("PkBlur") : 0;
  _smearPc = store.find_dstore("PcBlur") ?
    store.fetch_dstore("PcBlur") : 0;
  
  _UseTimeCut = store.find_istore("UseTimeCut") ?
    (store.fetch_istore("UseTimeCut")==1 ? true : false) : false;
  
  _edges = store.fetch_vstore("edges");

  if(_doVary != 0) {
    std::string xsecFile = store.fetch_sstore("xsecFile");
    _xsecF = new TFile( xsecFile.c_str(), "read");
    _xsecF->GetObject("xsec", _xsecErr);
  }

  TString cuts = "Fitted == 1";
  cuts += " && Charge.recQ == ";
  cuts += _beamCharge;
  cuts += " && abs(0.001/Momentum.recQ) < ";
  cuts += _maxOver * _maxE;
  cuts += " && TrajVertex[0] < ";
  cuts += _edges[0];
  cuts += " && fittedNode[0]/InMu[0] > ";
  cuts += _minNodes;

  _cuts = cuts;

}

/********************************************/
void Run_TMVA::define_variables(const bhep::gstore& store)
/********************************************/
{
  // Now add variables to the reader
  // In principle this can be controlled by the configuration file as well.
  vector<int> selvars;
  if(store.find_ivstore("MVAVars"))
    selvars = store.fetch_ivstore("MVAVars");
  else { 
    for(int i=0; i<9; i++) selvars.push_back(1);
    selvars[3] = 0;
    selvars[6] = 0;
  }
  _reader->AddSpectator( "nuEnergy  := NeuEng * 0.001", &nuEnergy);
  if(selvars[0]) _reader->AddVariable( "ErrqP  := ErrMom[0]/RecMom[0]", &ErrqP);
  if(selvars[1]) _reader->AddVariable( "trHits := InMu[0]", &trHits);
  if(selvars[2]) _reader->AddVariable( "Rp := initrangqP[0]/RecMom[0]", &Rp);
  if(selvars[3]) _reader->AddVariable( "engfrac := visEngTraj[0]/visibleEng", &engfrac);
  if(selvars[4]) _reader->AddVariable( "meanDep := trajEngDep[0]", &meanDep);
  if(selvars[5]) _reader->AddVariable( "EngVar := trajLowDep[0]/trajHighDep[0]", &EngVar);   
  if(selvars[6]) _reader->AddVariable( "recMom    := abs( 0.001/RecMom[0] )", &recMom);
  // if(selvars[7]) _reader->AddVariable( "Qt  :=  abs( 0.001/RecMom[0] ) * ( 1 - pow( hadRdirP[0][0]*RecXDirection[0] + hadRdirP[0][1]*RecYDirection[0] + hadRdirP[0][2], 2)/(1 + pow(RecXDirection[0],2) + pow(RecYDirection[0],2)))", &Qt );
  if(selvars[7]) _reader->AddVariable( "Qt  :=  abs( 0.001/RecMom[0] ) * ( 1 - pow(RecShowerDir[0]*RecXDirection[0] + RecShowerDir[1]*RecYDirection[0] + RecShowerDir[2], 2)/(1 + pow(RecXDirection[0],2) + pow(RecYDirection[0],2)))", &Qt );
  if(selvars[8]) _reader->AddVariable( "meanHPP := Nallhits/Planes", &meanHPP);


  TString dir    = "weights/";
  TString prefix = store.find_sstore("TrainingLib") ?
    store.fetch_sstore("TrainingLib") : "NC_Classification";

  // Book method(s)
  for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
    if (it->second) {
      TString methodName = TString(it->first) + TString(" method");
      TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
      _reader->BookMVA( methodName, weightfile );
    }
  }
}

void Run_TMVA::define_histograms(const bhep::gstore& store) {

  _OutFile->cd();
  int ntBins = store.fetch_istore("truBins");
  double trEdge[ntBins+1], recEdge[ntBins+2];
  define_rec_binning( store, recEdge, trEdge );

  alltracks = new TH1D("alltracks",";True Neutrino Energy (GeV)", ntBins, trEdge);
  alltracks2 = new TH2D("alltracks2",
			";Reconstructed Neutrino Energy (GeV);True Neutrino Energy (GeV)", 
			ntBins+1, recEdge, ntBins, trEdge);
  UInt_t nbin = 100; 
  pcuts = new TH1D("cuts_passed",";Events Passed by Cuts",11,0.0,11.0);
  rcuts = new TH1D("cuts_reject",";Events Rejected by Cuts",10,0.0,10.0);
  ntrajectories = new TH1D("ntrajectories",";Number of Trajectories Found",10,0,10); 

  CreateHistogram("No_MVA",        ntBins, recEdge, trEdge, nbin, -1, 1 );
  if (Use["Likelihood"])    CreateHistogram("Likelihood",    ntBins, recEdge, trEdge, nbin, -1, 1 );
  if (Use["LikelihoodD"])   CreateHistogram("LikelihoodD",   ntBins, recEdge, trEdge, nbin, -1, 0.9999 );
  if (Use["LikelihoodPCA"]) CreateHistogram("LikelihoodPCA", ntBins, recEdge, trEdge, nbin, -1, 1 );
  if (Use["LikelihoodKDE"]) CreateHistogram("LikelihoodKDE", ntBins, recEdge, trEdge, nbin,  -0.00001, 0.99999 );
  if (Use["LikelihoodMIX"]) CreateHistogram("LikelihoodMIX", ntBins, recEdge, trEdge, nbin,  0, 1 );
  if (Use["PDERS"])         CreateHistogram("PDERS",         ntBins, recEdge, trEdge, nbin,  0, 1 );
  if (Use["PDERSD"])        CreateHistogram("PDERSD",        ntBins, recEdge, trEdge, nbin,  0, 1 );
  if (Use["PDERSPCA"])      CreateHistogram("PDERSPCA",      ntBins, recEdge, trEdge, nbin,  0, 1 );
  if (Use["KNN"])           CreateHistogram("KNN",           ntBins, recEdge, trEdge, nbin, 0, 1 );
  if (Use["HMatrix"])       CreateHistogram("HMatrix",       ntBins, recEdge, trEdge, nbin, -0.95, 1.55 );
  if (Use["Fisher"])        CreateHistogram("Fisher",        ntBins, recEdge, trEdge, nbin, -4, 4 );
  if (Use["FisherG"])       CreateHistogram("FisherG",       ntBins, recEdge, trEdge, nbin, -1, 1 );
  if (Use["BoostedFisher"]) CreateHistogram("BoostedFisher", ntBins, recEdge, trEdge, nbin, -2, 2 );
  if (Use["LD"])            CreateHistogram("LD",            ntBins, recEdge, trEdge, nbin, -2, 2 );
  if (Use["MLP"])           CreateHistogram("MLP",           ntBins, recEdge, trEdge, nbin, -1.25, 1.5 );
  if (Use["MLPBFGS"])       CreateHistogram("MLPBFGS",       ntBins, recEdge, trEdge, nbin, -1.25, 1.5 );
  if (Use["MLPBNN"])        CreateHistogram("MLPBNN",        ntBins, recEdge, trEdge, nbin, -1.25, 1.5 );
  if (Use["CFMlpANN"])      CreateHistogram("CFMlpANN",      ntBins, recEdge, trEdge, nbin, 0, 1 );
  if (Use["TMlpANN"])       CreateHistogram("TMlpANN",       ntBins, recEdge, trEdge, nbin, -1.3, 1.3 );
  if (Use["BDT"])           CreateHistogram("BDT",           ntBins, recEdge, trEdge, nbin, -0.8, 0.8 );
  if (Use["BDTD"])          CreateHistogram("BDTD",          ntBins, recEdge, trEdge, nbin, -0.8, 0.8 );
  if (Use["BDTG"])          CreateHistogram("BDTG",          ntBins, recEdge, trEdge, nbin, -1.0, 1.0 );
  if (Use["RuleFit"])       CreateHistogram("RuleFit",       ntBins, recEdge, trEdge, nbin, -2.0, 2.0 );
  if (Use["SVM_Gauss"])     CreateHistogram("SVM_Gauss",     ntBins, recEdge, trEdge, nbin,  0.0, 1.0 );
  if (Use["SVM_Poly"])      CreateHistogram("SVM_Poly",      ntBins, recEdge, trEdge, nbin,  0.0, 1.0 );
  if (Use["SVM_Lin"])       CreateHistogram("SVM_Lin",       ntBins, recEdge, trEdge, nbin,  0.0, 1.0 );
  if (Use["FDA_MT"])        CreateHistogram("FDA_MT",        ntBins, recEdge, trEdge, nbin, -2.0, 3.0 );
  if (Use["FDA_GA"])        CreateHistogram("FDA_GA",        ntBins, recEdge, trEdge, nbin, -2.0, 3.0 );
  if (Use["Category"])      CreateHistogram("Category",      ntBins, recEdge, trEdge, nbin, -2., 2. );
  if (Use["Plugin"])        CreateHistogram("Plugin",        ntBins, recEdge, trEdge, nbin, -0.8, 0.8 );
  /* 
  // PDEFoam also returns per-event error, fill in histogram, and also fill significance
  if (Use["PDEFoam"]) {
    histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
    histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
    histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
  }

  // Book example histogram for probability (the other methods are done similarly)
  TH1F *probHistFi(0), *rarityHistFi(0);
  if (Use["Fisher"]) {
    probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
    rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
  }
  */
}
void Run_TMVA::CreateHistogram(std::string Key, int ntBins, double* recEdge, double* trEdge, 
		     int nbins, double min, double max){
  std::string histname = "MVA_"+Key;
  double cbins[nbins];
  for(int i=0;i<=nbins;i++) cbins[i] = min + double(i)*(max - min)/double(nbins);
  hcut[Key] = new TH1F( histname.c_str(), histname.c_str(), nbins, min, max);
  histname = "Eff"+Key;
  hEff[Key] = new TH1D( histname.c_str(),";True Energy (GeV)",ntBins, trEdge);
  histname = "recEff"+Key;
  hrecEff[Key] = new TH2D( histname.c_str(),
			   ";Reconstructed Energy (GeV);True Energy (GeV)",
			   ntBins+1, recEdge, ntBins, trEdge);
  histname = "recEff3D_"+Key;
  hrecEff3[Key] = new TH3D( histname.c_str(),
			    ";Reconstructed Energy (GeV);True Energy (GeV); MVA Cut",
			    ntBins+1, recEdge, ntBins, trEdge, nbins, cbins);
  
  histname = "Nhits_" + Key;
  hNhits[Key] = new TH1D(histname.c_str(),";Number of hits in Trajectory", 140, 0, 280);
  histname = "ErrqP_" + Key;
  hErrqP[Key] = new TH1D(histname.c_str(),";#sigma_{q/p}/(q/p)", 1000, -5, 5);
  histname = "Rp_" + Key;
  hRp[Key]    = new TH1D(histname.c_str(),";(q_{init}/p_{range})#times(p_{fit}/q_{fit})",100,-4,4);
  histname = "MeanEDep_" + Key;
  hEngDep[Key]= new TH1D(histname.c_str(),";#sum #Delta E_{low} / #sum #Delta E_{high}",450, 0, 45.);
  histname = "EngVar_" + Key;
  hEngVar[Key]= new TH1D(histname.c_str(),";Mean #Delta E (MeV)",100, 0, 1);
  histname = "Qt_" + Key;
  hQt[Key]    = new TH1D(histname.c_str(),";Energy Transfer Q_{t} (GeV)",100, 0, 10);
  histname = "meanHPP_" + Key;
  hmeanHPP[Key] = new TH1D(histname.c_str(),";Mean Hits Per Plane", 100, 0, 20);
  histname = "recEff_QES_" + Key;
  hrecEffQES[Key] = new TH2D( histname.c_str(),
			      ";Reconstructed Energy (GeV);True Energy (GeV)",
			      ntBins+1, recEdge, ntBins, trEdge);
  histname = "recEff_DIS_" + Key;
  hrecEffDIS[Key] = new TH2D( histname.c_str(),
			      ";Reconstructed Energy (GeV);True Energy (GeV)",
			      ntBins+1, recEdge, ntBins, trEdge);
  histname = "recEff_nQES_" + Key;
  hrecEffnQES[Key] = new TH2D( histname.c_str(),
			      ";Reconstructed Energy (GeV);True Energy (GeV)",
			      ntBins+1, recEdge, ntBins, trEdge);
  histname = "recEff_nDIS_" + Key;
  hrecEffnDIS[Key] = new TH2D( histname.c_str(),
			      ";Reconstructed Energy (GeV);True Energy (GeV)",
			      ntBins+1, recEdge, ntBins, trEdge);
  histname = "timeRev_" + Key;
  htimeRev[Key] = new TH2D( histname.c_str(),
  			    ";Reconstructed Momentum (GeV); True Momentum (GeV)",
  			    100, -10, 10, 100, -10, 10);

}

void Run_TMVA::FillHistogram(std::string Key, double Erec, double Etru,
			     double qPrec, double qPtru, 
			     double cutvar, double limit, int trutype){
  hcut[Key]->Fill(cutvar);
  int q = qPrec != 0.0 ? int(qPrec/fabs(qPrec)) : 0;
  if(Erec > _maxE) Erec = _maxE;
  if(q == _beamCharge){
    if(cutvar > limit){
      hEff[Key]->Fill(Etru);
      hrecEff[Key]->Fill(Erec, Etru);
      if(trutype==1 || trutype==2) hrecEffQES[Key]->Fill(Erec, Etru);
      if(trutype==3 || trutype==4) hrecEffDIS[Key]->Fill(Erec, Etru);
      if(trutype!=1 && trutype!=2) hrecEffnQES[Key]->Fill(Erec, Etru);
      if(trutype!=3 && trutype!=4) hrecEffnDIS[Key]->Fill(Erec, Etru);
      hNhits[Key]->Fill(trHits);
      hRp[Key]->Fill(Rp);
      hErrqP[Key]->Fill(ErrqP);
      hEngDep[Key]->Fill(meanDep);
      hEngVar[Key]->Fill(EngVar);
      hQt[Key]->Fill(Qt);
      hmeanHPP[Key]->Fill(meanHPP);
      htimeRev[Key]->Fill(0.001/qPrec, 0.001/qPtru);
    }
    hrecEff3[Key]->Fill(Erec, Etru, cutvar); 
  }
}
/********************************************/
void Run_TMVA::plotEventTopology(int ifig, int imu, std::string method, 
				 std::vector<std::vector<double> >* xpos,
				 std::vector<std::vector<double> >* ypos,
				 std::vector<std::vector<double> >* zpos){
  if( xpos->at(imu).size() == ypos->at(imu).size() && 
      xpos->at(imu).size() == zpos->at(imu).size()){
    geventDisplayRZ[method].push_back(new TGraph(xpos->at(imu).size()));
    std::ostringstream s1;
    s1<<"Topology_"<<method<<"_RZ_"<<ifig;
    std::string rzname(s1.str());
    geventDisplayRZ[method].back()->SetNameTitle(rzname.c_str(), 
						 ";Z Position of Hit (metres); Radial Position of Hit (metres)"); 
    geventDisplayYX[method].push_back(new TGraph(xpos->at(imu).size()));
    
    std::ostringstream s2;
    s2<<"Topology_"<<method<<"_YX_"<<ifig;
    std::string yxname(s2.str());
    geventDisplayYX[method].back()->SetNameTitle(yxname.c_str(), 
						 ";X Position of Hit (metres); Y Position of Hit (metres)"); 
    for( int ip=0; ip<xpos->at(imu).size(); ip++){
      double R = sqrt(xpos->at(imu).at(ip)*xpos->at(imu).at(ip) + 
		      ypos->at(imu).at(ip)*ypos->at(imu).at(ip))/1000.;
      geventDisplayRZ[method].back()->SetPoint(ip,zpos->at(imu).at(ip)/1000.,R);
      geventDisplayYX[method].back()->SetPoint(ip,xpos->at(imu).at(ip)/1000.,
					       ypos->at(imu).at(ip)/1000.);
    }
  }
}


/********************************************/
void Run_TMVA::define_rec_binning(const bhep::gstore& store, double* recB, double* truB)
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
void Run_TMVA::execute()
/********************************************/
{
  
  TStopwatch sw;
  sw.Start();

  // Do the analysis
  theTree = new TChain("tree");
  Int_t planes, nallhits, imu=0;
  Int_t   fit[10], no_traj, truInt;
  Int_t inmu[10], fitn[10], Q[10], hadhits;
  Double_t qP[10], errqP[10], initqP[10], NeuEng, visibleEng, qPtru[10];
  Double_t visEngTraj[10], trajEngDep[10], trajLowDep[10], trajHighDep[10];
  Double_t QESEng[2], truHad, hadMom[3], hadRdirP[2][3], muXdir[10], muYdir[10];
  Double_t nonmuEDep, recShowerDir[3];
  Double_t vertX[10], vertY[10], vertZ[10], truvert[3];

  TBranch *b_XPositions;
  TBranch *b_YPositions;
  TBranch *b_ZPositions;
  TBranch *b_HTime;

  _XPos = 0;
  _YPos = 0;
  _ZPos = 0;
  _HTime = 0;

  for ( int ifile = _firstT; ifile <= _lastT; ifile++){
    std::string treefileName = _filebaseTree+bhep::to_string( ifile )+".root";
    std::cout<<treefileName<<std::endl;
    TFile* treeFile = new TFile( treefileName.c_str() );
    if(treeFile->IsZombie()) continue;
    
    treeFile->Close();
  
    // TTree* theTree;
    // treeFile->GetObject("tree", theTree);
    theTree->Add(treefileName.c_str());
  }
  // theTree->Print();
  theTree->SetBranchStatus("*",0);
  
  theTree->SetBranchStatus( "TrajectoryNo", 1);    theTree->SetBranchAddress("TrajectoryNo", &no_traj);
  theTree->SetBranchStatus( "Fitted", 1 );         theTree->SetBranchAddress( "Fitted", &fit );
  theTree->SetBranchStatus( "TrueInteraction", 1); theTree->SetBranchAddress( "TrueInteraction", &truInt );
  theTree->SetBranchStatus( "TruMom", 1 );         theTree->SetBranchAddress( "TruMom", &qPtru );
  theTree->SetBranchStatus( "RecMom", 1 );         theTree->SetBranchAddress( "RecMom", &qP );
  theTree->SetBranchStatus( "ErrMom", 1 );         theTree->SetBranchAddress( "ErrMom", &errqP );
  theTree->SetBranchStatus( "RecCharge", 1 );      theTree->SetBranchAddress( "RecCharge", &Q );
  // theTree->SetBranchStatus( "LongMuTraj", 1 );     theTree->SetBranchAddress("LongMuTraj", &imu );
  theTree->SetBranchStatus( "NeuEng", 1 );         theTree->SetBranchAddress( "NeuEng", &NeuEng );
  theTree->SetBranchStatus( "Planes", 1);          theTree->SetBranchAddress( "Planes", &planes );
  theTree->SetBranchStatus( "Nallhits",1);         theTree->SetBranchAddress( "Nallhits", &nallhits);
  theTree->SetBranchStatus( "InMu", 1 );           theTree->SetBranchAddress( "InMu", &inmu );
  if( theTree->GetBranch("FitN") ){
    theTree->SetBranchStatus( "FitN", 1 );     theTree->SetBranchAddress( "FitN", &fitn );
  } else {
    theTree->SetBranchStatus( "fittedNode", 1 );     theTree->SetBranchAddress( "fittedNode", &fitn );
  }
  theTree->SetBranchStatus( "initrangqP", 1);      theTree->SetBranchAddress( "initrangqP", &initqP);
  theTree->SetBranchStatus( "visibleEng", 1 );     theTree->SetBranchAddress( "visibleEng", &visibleEng );
  theTree->SetBranchStatus( "visEngTraj", 1 );     theTree->SetBranchAddress( "visEngTraj", &visEngTraj );
  theTree->SetBranchStatus( "trajEngDep", 1 );     theTree->SetBranchAddress( "trajEngDep", &trajEngDep );
  theTree->SetBranchStatus( "trajLowDep", 1 );     theTree->SetBranchAddress( "trajLowDep", &trajLowDep );
  theTree->SetBranchStatus( "trajHighDep",1 );     theTree->SetBranchAddress( "trajHighDep",&trajHighDep);
  // theTree->SetBranchStatus( "hadQESEng", 1 );      theTree->SetBranchAddress( "hadQESEng", &QESEng );
  // theTree->SetBranchStatus( "hadRecEng", 1 );      theTree->SetBranchAddress( "hadRecEng", &hadDep );
  theTree->SetBranchStatus( "hadTrueE", 1 );      theTree->SetBranchAddress( "hadTrueE", &truHad );
  theTree->SetBranchStatus( "hadronP", 1 );       theTree->SetBranchAddress( "hadronP", &hadMom );
  theTree->SetBranchStatus( "hadRdirP",  1 );      theTree->SetBranchAddress( "hadRdirP",  &hadRdirP);
  theTree->SetBranchStatus( "NonMuonEdep", 1 );      theTree->SetBranchAddress( "NonMuonEdep", &nonmuEDep );
  theTree->SetBranchStatus( "hadronHits", 1 );      theTree->SetBranchAddress( "hadronHits", &hadhits );
  theTree->SetBranchStatus( "RecShowerDir",  1 );      theTree->SetBranchAddress( "RecShowerDir",  &recShowerDir);
  theTree->SetBranchStatus( "RecXDirection", 1);   theTree->SetBranchAddress( "RecXDirection", &muXdir);
  theTree->SetBranchStatus( "RecYDirection", 1);   theTree->SetBranchAddress( "RecYDirection", &muYdir);
  theTree->SetBranchStatus("TrajVertex",1);       theTree->SetBranchAddress("TrajVertex", &vertZ);   
  theTree->SetBranchStatus("evVertex",1);       theTree->SetBranchAddress("evVertex", &truvert);   
  theTree->SetBranchStatus("RecXPosition",1);    theTree->SetBranchAddress("RecXPosition", &vertX);
  theTree->SetBranchStatus("RecYPosition",1);    theTree->SetBranchAddress("RecYPosition", &vertY);
  // theTree->SetBranchStatus( "Direction", 1 );theTree->SetBranchAddress( "Direction", &Dir );
  if( theTree->GetBranch( "XPositions" ) ){
    theTree->SetBranchStatus("XPositions",1);       theTree->SetBranchAddress("XPositions", &_XPos, &b_XPositions);}
  else{
    theTree->SetBranchStatus("XPos",1);       theTree->SetBranchAddress("XPos", &_XPos, &b_XPositions);}
  if( theTree->GetBranch( "YPositions" ) ){
    theTree->SetBranchStatus("YPositions",1);       theTree->SetBranchAddress("YPositions", &_YPos, &b_YPositions);}
  else{
    theTree->SetBranchStatus("YPos",1);       theTree->SetBranchAddress("YPos", &_YPos, &b_YPositions);}
  if( theTree->GetBranch( "ZPositions" ) ){
    theTree->SetBranchStatus("ZPositions",1);       theTree->SetBranchAddress("ZPositions", &_ZPos, &b_ZPositions);}
  else{
    theTree->SetBranchStatus("ZPos",1);       theTree->SetBranchAddress("ZPos", &_ZPos, &b_ZPositions);}
  if( theTree->GetBranch("Time") ) {
    theTree->SetBranchStatus("Time",1);       theTree->SetBranchAddress("Time", &_HTime, &b_HTime);}
    
  // Efficiency calculator for cut method
  Int_t    nSelCutsGA = 0;
  Double_t effS       = 0.7;
  bool plotEvent      = false; 
  int  eventsPlotted  = 0;
  
  std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
  
  std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
  for (Long64_t ievt=0; ievt<theTree->GetEntries()-1;ievt++) {
    
    if (ievt%1000 == 0) 
      std::cout << "--- ... Processing event: " << ievt << std::endl;
    
    bool plotEvent=false;
    if(rand->Rndm() < 0.05){ plotEvent=true; }
    theTree->GetEntry(ievt);
    // std::cout<<"For Event "<<ievt<<" with Fitted = "<<fit<<", NeuEng = "<<NeuEng<<", RecMom[LongMuTraj] = "<<1./qP[1]<<std::endl;
    /*
    if(_doVary > 0){
      int bin = _xsecErr->FindBin( NeuEng/1000 );
      double xsErr = _xsecErr->GetBinContent( bin );
      if( truInt == _Vtype && _doVary == 1)
	// We are reducing the effective cross-section by 10%
	// Throw a random number
	if(rand->Rndm() <= 1-xsErr) 
	  continue;
      if( truInt != _Vtype && _doVary == 2)
	// We are reducing all cross sections by 10% except for target interaction
	// Throw a random number 
	if(rand->Rndm() < 1/(xsErr+1)) 
	  continue;
    }
    */
    bool PassQualityCuts = true;
    alltracks->Fill(NeuEng * 0.001);
    fill_tracks2();
    ntrajectories->Fill(no_traj);
    pcuts->Fill("Initial Events",1);
    if(no_traj == 0){
      rcuts->Fill("Trajectories Found",1);
      continue;
    } else pcuts->Fill("Trajectories Found",1);
    if(fit[imu]==0){
      rcuts->Fill("Long Trajectory Fit Success",1);
      continue;
    } else pcuts->Fill("Long Trajectory Fit Success",1);
    
    //if ( !fiducial(vertX[imu], vertY[imu], vertZ[imu]) ){
    if ( !fiducial2(no_traj)){
      rcuts->Fill("Fiducial",1);
      continue;
    } else pcuts->Fill("Fiducial",1); 
    /*
    if ( _UseTimeCut && !pseudoTimeCut(_XPos->at(imu), _YPos->at(imu), 
				      _ZPos->at(imu), qP[imu], qPtru[0])){
      rcuts->Fill("Pseudo-Time Cut",1);
      continue;
    } else if( _UseTimeCut )
      pcuts->Fill("Pseudo-Time Cut", 1);
    */

    nuEnergy = NeuEng * 0.001;
    double sigmaqP = _smearPr * (_smearPc + _smearPk*fabs(qP[imu]) + _smearPp/fabs(qP[imu])) * qP[imu];
    if(_smearPr > 0) qP[imu] +=  rand->Gaus(sigmaqP);
    ErrqP = errqP[imu]/qP[imu];
    trHits = inmu[imu];
    
    Rp = initqP[imu]/qP[imu];
    if(Rp < 0) {
      rcuts->Fill("Range Fit Ratio",1);
      continue;
    } else pcuts->Fill("Range Fit Ratio",1);
    
    engfrac = visEngTraj[imu]/visibleEng;
    meanHPP = float(nallhits)/float(planes);
    meanDep = trajEngDep[imu];
    EngVar = trajLowDep[imu]/trajHighDep[imu];
    recMom = fabs(0.001/qP[imu]);

    double recEng = recMom;
    double dxdz = muXdir[imu];
    double dydz = muYdir[imu];
    double mP  = 938.27 * MeV;
    double mmu = 105.65 * MeV;
    double mN  = 939.56 * MeV;
    double Emu = sqrt(1./qP[imu]/qP[imu] + mmu*mmu);
    double costhmu = 1./sqrt(1. + dxdz*dxdz + dydz*dydz);
    if(qP[imu] > 0) // assume the QE interaction numubar + p \to n + mu^{+}                                                                                                                       
      QESEng[0] = (mP * Emu + 0.5*(mN*mN - mmu*mmu - mP*mP))/(mP - Emu + costhmu/fabs(qP[imu]));
    else if(qP[imu] < 0) // assume the QE interaction numu + n \to p + mu^{-}                                                                                                                     
      QESEng[0] = (mN * Emu + 0.5*(mP*mP - mmu*mmu - mN*mN))/(mN - Emu + costhmu/fabs(qP[imu]));
    else
      QESEng[0] = 0.;

    
    if (  hadhits==0 ) 
      recEng = fabs(QESEng[0])*0.001;
    else if ( hadhits > 0 ){
      double altEng = NeuEng - sqrt(1./qPtru[0]/qPtru[0] + mmu*mmu);
      if(altEng > 0) altEng *= 0.001;
      else altEng = 0.00;
	// 0.001 * sqrt(pow(hadMom[0],2) + pow(hadMom[1],2) + pow(hadMom[2],2));
      // if ( _beamCharge != fabs(qP[imu])/qP[imu] )
      // altEng -= recMom;
      double Esig = sqrt(pow(_smearRoot*sqrt( altEng ),2) + pow(_smearConst*altEng,2));
      recEng = sqrt( recMom*recMom + 0.10565*0.10565) + altEng + rand->Gaus(0,Esig); 
    }
    else {
      recEng = 0.;
    }
    /*
    if (hadhits == 0 && fabs(QESEng[0]) > 0) recEng = fabs(QESEng[0])*0.001;
    else if (hadhits > 0){
      double Ehad =nonmuEDep * (1 + 1.5 / 0.135 * 1.451 / 1.936 ) * 0.001; 
      recEng = sqrt( recMom * recMom + 0.10566 * 0.10566 ) + Ehad;
      if(_smearRoot > 0.0){
	double Esig = sqrt(pow(_smearRoot*sqrt(altEng),2) + pow(_smearConst*altEng,2));
	recEng += rand->Gaus(0,Esig); 
    }
    else {
      recEng = 0.;
    }
    */
    double dirnorm2 = (1 + pow(muXdir[imu],2) + pow(muYdir[imu],2));
    double cos2th = pow( recShowerDir[0]*muXdir[imu] +
			 recShowerDir[1]*muYdir[imu] +
			 recShowerDir[2], 2)/dirnorm2;
    double sin2th = 1 - cos2th; 
    Qt = recMom * sin2th; 
    // std::cout<<"Fill the normalization histogram"<<std::endl;
    
    // std::cout<<"If the event fails some data quality cuts then reject the event prior to MVA"<<std::endl;
    
    // if(fit==0) PassQualityCuts = false;
    if(fabs(_beamCharge * 0.001/qP[imu]) > _maxOver * _maxE){
      rcuts->Fill("Maximum Momentum",1);
      PassQualityCuts = false;
    } else pcuts->Fill("Maximum Momentum",1);
    if(double(fitn[imu])/double(inmu[imu]) < _minNodes){
      rcuts->Fill("Fitted Proportion",1);
      PassQualityCuts = false;
    } else pcuts->Fill("Fitted Proportion",1);
    if( !PassQualityCuts ){
      //std::cout<<"Event has failed cuts: Fitted = "<<fit<<", Momentum cut = "
      //	 << fabs(0.001/qP[1]) <<":"<<_maxOver * _maxE
      //	 <<", Proportion Cut = "<<double(hit[3])<<":"<<hit[1]<<std::endl;
      continue;
    }
    std::string method = "No_MVA";
    double weight = 1.0;
    double weightcut = 0.5;
    FillHistogram(method, recEng, nuEnergy, qP[imu], qPtru[0], weight, weightcut, truInt); 
    for(int i=0; i<_myMethodList.size(); i++){
      std::string method = _myMethodList[i];
      double weightcut = _MVAcut[i];
      double weight = _reader->EvaluateMVA( method+" method") ;
      FillHistogram(method, recEng, nuEnergy, qP[imu], qPtru[0], weight, weightcut, truInt); 
      if( eventsPlotted < 100 && weight > weightcut){
	if(plotEvent){
	  plotEventTopology(eventsPlotted, imu, method, _XPos, _YPos, _ZPos);
	  eventsPlotted++;
	}
      }
      
      if ( weight > weightcut ){
	pcuts->Fill(method.c_str(),1);
      } else {
	rcuts->Fill(method.c_str(),1);
      }
      
    }
  }
  
  
  std::cout<<" Get elapsed time"<<std::endl;
  sw.Stop();
  std::cout << "--- End of event loop: "; sw.Print();
  
}

/********************************************/
void Run_TMVA::fill_tracks2()
/********************************************/
{
  //Fill a column of the 2D matrix for normalisation.
  double incr = 0.0;
  int bin;
  while ( incr < alltracks2->GetXaxis()->GetXmax() )
    {
      alltracks2->Fill( incr, nuEnergy );

      bin = alltracks2->GetXaxis()->FindBin( incr );
      incr += alltracks2->GetXaxis()->GetBinWidth( bin );
    }

}
/********************************************/
bool Run_TMVA::fiducial(double vx, double vy, double vz){
/********************************************/
  if ( vz > _edges[0] ) return false;
  if ( vz < _edges[1] ) return false;
  double vertR = sqrt(vx*vx + vy*vy);
  if ( vertR > _edges[2] ) return false;
  if ( vertR < _edges[3] ) return false;
  if( fabs(vx) > _edges[4] ) return false;
  if( fabs(vy) > _edges[5] ) return false;
  if( fabs(vy) > _edges[4] * (1 + tan(atan(1.)/2.)) - fabs(vx) ) return false;
  return true;
}
/********************************************/
bool Run_TMVA::fiducial2(int nTraj) { 
/********************************************/
  // Start by finding the maximum and minimum z for the track
  double minx=999999, maxx=-999999;
  double miny=999999, maxy=-999999;
  double minz=999999, maxz=-999999;
  double mint=999999, maxt=-999999;
  double sumdt = 0;
  for(int itrack=0; itrack<nTraj; itrack++){
    std::vector<double>::iterator xit=_XPos->at(itrack).begin();
    std::vector<double>::iterator yit=_YPos->at(itrack).begin();
    std::vector<double>::iterator zit=_ZPos->at(itrack).begin();
    std::vector<double>::iterator timeit=_HTime->at(itrack).begin();
    std::vector<double>::iterator it;
    std::vector<double>::iterator term;
    //if(Time.size() > 0){
    //  it=Time.begin();
    //  term=Time.end();
    //}
    //else{
    it=_ZPos->at(itrack).begin();
    term=_ZPos->at(itrack).end();
    //}
    double lastz=0.0, currentz=0.0;
    double lastt=0.0, currentt=0.0;
    double sumvz=0.0;
    for(it;it!=term;it++){
      if(it!=_ZPos->at(0).begin() && (*timeit) > 0.05){
	lastz=currentz;
	lastt=currentt;
      }
      currentz=(*zit);
      currentt=(*timeit);
      if(it!=_ZPos->at(0).begin() && (*timeit) > 0.05)
	sumvz=(currentz - lastz)/(currentt-lastt);
      
      if((*it) < minz){
	minx = (*xit);
	miny = (*yit);
	minz = (*zit);
	if( (*timeit) > 0.05 ) mint = (*timeit);
      }
      if((*it) > maxz){
	maxx = (*xit);
	maxy = (*yit);
	maxz = (*zit);
	if( (*timeit) > 0.05 ) maxt = (*timeit);
      }
      xit++;
      yit++;
      zit++;
      timeit++;
    }
    
    double vertRmin = sqrt(minx*minx + miny*miny);
    double vertRmax = sqrt(maxx*maxx + maxy*maxy);
    if (sumvz < 0) { // maxt < mint){
      // switch the first and last points
      double sz=maxz, sy=maxy, sx=maxx, sr=vertRmax, st=maxt;
      maxx=minx;
      maxy=miny;
      maxz=minz;
      maxt=mint;
      vertRmax = vertRmin;
      minx=sx;
      miny=sy;
      minz=sz;
      mint=st;
      vertRmin = sr;
    }
    if (      minz > _edges[0] ) return false;
    if (      minz < _edges[1] ) return false;
    if (      maxz > _edges[0] ) return false;
    if (      maxz < _edges[1] ) return false;
    if (  vertRmin > _edges[2] ) return false;
    if (  vertRmin < _edges[3] ) return false;
    if((fabs(minx)) > _edges[4] ) return false;
    if((fabs(miny)) > _edges[5] ) return false;
    if((fabs(miny) + fabs(minx)) > _edges[4] * (1 + tan(atan(1.)/2.)) ) return false;
    if((fabs(maxx)) > _edges[4] ) return false;
    if((fabs(maxy)) > _edges[5] ) return false;
    if((fabs(maxy) + fabs(maxx)) > _edges[4] * (1 + tan(atan(1.)/2.)) ) return false;
    // std::cout<<"("<<minx<<", "<<miny<<", "<<minz<<", "<<vertRmin<<", "<<mint<<") to ("
    //	   <<maxx<<", "<<maxy<<", "<<maxz<<", "<<vertRmax<<", "<<maxt<<") with sumvz = "<<sumvz<<"\n";
  }
  return true;
  
}
/********************************************/
bool Run_TMVA::pseudoTimeCut(std::vector<double> XPos,
			     std::vector<double> YPos,
			     std::vector<double> ZPos,
			     double qP, double qPtru){
  // Start by finding the maximum and minimum z for the track
  double minx=999999, maxx=-999999;
  double miny=999999, maxy=-999999;
  double minz=999999, maxz=-999999;
  std::vector<double>::iterator xit=XPos.begin();
  std::vector<double>::iterator yit=YPos.begin();
  std::vector<double>::iterator zit;
  for(zit=ZPos.begin();zit!=ZPos.end();zit++){
    if((*zit) < minz){
      minx = (*xit);
      miny = (*yit);
      minz = (*zit);
    }
    if((*zit) > maxz){
      maxx = (*xit);
      maxy = (*yit);
      maxz = (*zit);
    }
  }
  double maxR = sqrt(maxx*maxx + maxy*maxy);
  double minR = sqrt(minx*minx + miny*miny);
  // Now check the relative charge of the reconstruction
  double chargeproduct = (qP * qPtru) / (fabs(qP) * fabs(qPtru));
  if( chargeproduct == 1){ // No discrepancy. No need to do more.
    if(minR > _edges[2]) return false;
    else if(minz > _edges[0]) return false;
    else return true;
  } else { // Check the position of the endpoint.
    //
    //
    // Is the maximum Z outside of the fiducial region
    /*
    if(maxR > _edges[2]) return false;
    else if(minR > _edges[2]) return false;
    else if(maxz > _edges[0]) return false;
    else if(minz < _edges[1]) return false;
    else return true;
    */
    // The above cut suggests that a track moving backwards is acceptable 
    // and requires that the "event vertex" is inside the fiducial. In practice
    // we know that all cosmic rays are inacceptable and can be cut.
    // That said this cut should never be used for neutrino simulations.
    return false;
  }
}

/********************************************/
void Run_TMVA::finalize()
/********************************************/
{
  // --- Write histograms
  _OutFile->cd();
  // _OutFile->Write();
 
  alltracks->Write();
  alltracks2->Write();
  rcuts->Write();
  pcuts->Write();
  ntrajectories->Write();
  for (int i=0; i<_myMethodList.size(); i++){
    hcut[_myMethodList[i]]->Write();
    hEff[_myMethodList[i]]->Sumw2();
    hEff[_myMethodList[i]]->Divide(alltracks);
    hEff[_myMethodList[i]]->Write();
    hrecEff[_myMethodList[i]]->Sumw2();
    hrecEff[_myMethodList[i]]->Divide(alltracks2);
    hrecEff[_myMethodList[i]]->Write();
    hrecEff3[_myMethodList[i]]->Write();
    hNhits[_myMethodList[i]]->Write();
    hRp[_myMethodList[i]]->Write();
    hErrqP[_myMethodList[i]]->Write();
    hEngDep[_myMethodList[i]]->Write();
    hEngVar[_myMethodList[i]]->Write();
    hQt[_myMethodList[i]]->Write();
    hmeanHPP[_myMethodList[i]]->Write();
    hrecEffQES[_myMethodList[i]]->Sumw2();
    hrecEffQES[_myMethodList[i]]->Divide(alltracks2);
    hrecEffQES[_myMethodList[i]]->Write();
    hrecEffDIS[_myMethodList[i]]->Sumw2();
    hrecEffDIS[_myMethodList[i]]->Divide(alltracks2);
    hrecEffDIS[_myMethodList[i]]->Write();
    hrecEffnQES[_myMethodList[i]]->Sumw2();
    hrecEffnQES[_myMethodList[i]]->Divide(alltracks2);
    hrecEffnQES[_myMethodList[i]]->Write();
    hrecEffnDIS[_myMethodList[i]]->Sumw2();
    hrecEffnDIS[_myMethodList[i]]->Divide(alltracks2);
    hrecEffnDIS[_myMethodList[i]]->Write();
    htimeRev[_myMethodList[i]]->Write();
    for(int iplot=0; iplot<geventDisplayRZ[_myMethodList[i]].size();  iplot++){
      geventDisplayRZ[_myMethodList[i]][iplot]->Write();
      geventDisplayYX[_myMethodList[i]][iplot]->Write();
    }
  }
  
  std::string _method = "No_MVA";
  
  hcut[_method]->Write();
  hEff[_method]->Sumw2();
  hEff[_method]->Divide(alltracks);
  hEff[_method]->Write();
  hrecEff[_method]->Sumw2();
  hrecEff[_method]->Divide(alltracks2);
  hrecEff[_method]->Write();
  hrecEff3[_method]->Write();
  hNhits[_method]->Write();
  hRp[_method]->Write();
  hErrqP[_method]->Write();
  hEngDep[_method]->Write();
  hEngVar[_method]->Write();
  hQt[_method]->Write();
  hmeanHPP[_method]->Write();
  hrecEffQES[_method]->Sumw2();
  hrecEffQES[_method]->Divide(alltracks2);
  hrecEffQES[_method]->Write();
  hrecEffDIS[_method]->Sumw2();
  hrecEffDIS[_method]->Divide(alltracks2);
  hrecEffDIS[_method]->Write();
  hrecEffnQES[_method]->Sumw2();
  hrecEffnQES[_method]->Divide(alltracks2);
  hrecEffnQES[_method]->Write();
  hrecEffnDIS[_method]->Sumw2();
  hrecEffnDIS[_method]->Divide(alltracks2);
  hrecEffnDIS[_method]->Write();
  htimeRev[_method]->Write();

  _OutFile->Close();
  delete _OutFile;
  delete _reader;
}

  
  /*
  if (Use["Likelihood"   ])   histLk     ->Write();
  if (Use["LikelihoodD"  ])   histLkD    ->Write();
  if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
  if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
  if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
  if (Use["PDERS"        ])   histPD     ->Write();
  if (Use["PDERSD"       ])   histPDD    ->Write();
  if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
  if (Use["KNN"          ]){
    histKNN    ->Write();
    recEffKNN  ->Sumw2();
    recEffKNN  ->Divide(alltracks2);
    recEffKNN  ->Write();
    recBackKNN ->Write();
    recEffKNNe ->Write();
    recBackKNNe->Write();
  }
  if (Use["HMatrix"      ])   histHm     ->Write();
  if (Use["Fisher"       ])   histFi     ->Write();
  if (Use["FisherG"      ])   histFiG    ->Write();
  if (Use["BoostedFisher"])   histFiB    ->Write();
  if (Use["MLP"           ]){
    histNn     ->Write();
    recEffNn   ->Sumw2();
    recEffNn   ->Divide(alltracks2);
    recEffNn   ->Write();
    recBackNn  ->Write();
    recEffNne  ->Write();
    recBackNne ->Write();
  }
  if (Use["LD"           ])   histLD     ->Write();
  if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
  if (Use["MLPBNN"       ])   histNnbnn  ->Write();
  if (Use["CFMlpANN"     ])   histNnC    ->Write();
  if (Use["TMlpANN"      ])   histNnT    ->Write();
  if (Use["BDT"          ]){
    histBdt    ->Write();
    recEffBdt  ->Sumw2();
    recEffBdt  ->Divide(alltracks2);
    recEffBdt  ->Write();
    recBackBdt ->Write();
    recEffBdte ->Write();
    recBackBdte->Write();
  }
  if (Use["BDTD"         ])   histBdtD   ->Write();
  if (Use["BDTG"         ])   histBdtG   ->Write(); 
  if (Use["RuleFit"      ])   histRf     ->Write();
  if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
  if (Use["SVM_Poly"     ])   histSVMP   ->Write();
  if (Use["SVM_Lin"      ])   histSVML   ->Write();
  if (Use["FDA_MT"       ])   histFDAMT  ->Write();
  if (Use["FDA_GA"       ])   histFDAGA  ->Write();
  if (Use["Category"     ])   histCat    ->Write();
  if (Use["Plugin"       ])   histPBdt   ->Write();
  
  // Write also error and significance histos
  if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }
  
  // Write also probability hists
  // if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
  
  
  */


//Test for voxel creation.

#include <mind/hit_clusterer.h>
#include <mind/fitter.h>

#include <bhep/gstore.h>
#include <bhep/hit.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>

#include <recpack/Measurement.h>
#include <recpack/stc_tools.h>

#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>

using namespace bhep;

int main(){

  int nhit1, nhit2, evt, nVox[2500], nhit3;
  double x1[5000], y1[5000], z1[5000];
  double xR[7000], yR[7000], zR[7000];
  double x2[5000], y2[5000], z2[5000], xRes[5000], yRes[5000], vE[5000], xVe[2][5000], rawE[5000],nvoxX[5000], nvoxY[5000];
  TFile *f1 = new TFile("mubarPlots.root","recreate");
  TTree *tree = new TTree("h1","stuff");
  tree->Branch("Event", &evt,"evtno/I");
  tree->Branch("G4Hit", &nhit3,"nhit3/I");
  tree->Branch("G4X",&xR,"xR[nhit3]/D");
  tree->Branch("G4Y",&yR,"yR[nhit3]/D");
  tree->Branch("G4Z",&zR,"zR[nhit3]/D");
  tree->Branch("bhepHit",&nhit1,"nhit1/I");
  tree->Branch("recHit",&nhit2,"nhit2/I");
  tree->Branch("bhepX",&x1,"x1[nhit1]/D");
  tree->Branch("bhepY",&y1,"y1[nhit1]/D");
  tree->Branch("bhepZ",&z1,"z1[nhit1]/D");
  tree->Branch("voxE", &vE,"vE[nhit1]/D");
  tree->Branch("endXE",&xVe[0],"xVE[nhit1]/D");
  tree->Branch("endYE",&xVe[1],"yVE[nhit1]/D");
  tree->Branch("RawVoxE", &rawE,"rawE[nhit1]/D");
  tree->Branch("recX",&x2,"x2[nhit2]/D");
  tree->Branch("recY",&y2,"y2[nhit2]/D");
  tree->Branch("recZ",&z2,"z2[nhit2]/D");
  tree->Branch("xResolution",&xRes,"xRes[nhit2]/D");
  tree->Branch("yResolution",&yRes,"yRes[nhit2]/D");
  tree->Branch("nVox",&nVox,"nvox[nhit2]/I");
  tree->Branch("nVoxX",&nvoxX,"nvoxX[nhit2]/D");
  tree->Branch("nVoxY",&nvoxY,"nvoxY[nhit2]/D");
  TH1F *voxXres = new TH1F("resX","Difference in voxel X and contained deposit Xs",200,-20,20);
  TH1F *voxYres = new TH1F("resY","Difference in voxel Y and contained deposit Ys",200,-20,20);
  int nevents=40000;

  string param_file = "examples/param/mubarCC_0.param";

  bhep::gstore store;
  bhep::sreader reader(store);
  reader.file(param_file);
  reader.group("RUN");
  reader.read();

  hit_clusterer* hct = new hit_clusterer( store );

  bhep::reader_root inDst;

  inDst.open("../digi_out/muCC/muCC_0_digi.dst.root");

  //measurement_vector meas;
  std::vector<cluster*> meas;
  size_t countification;
  for (int i=0;i<nevents;i++){
    cout << "Event: " << i << endl;
    
    bhep::event& e = inDst.read_event( i );
    
    std::vector<bhep::particle*> parts = e.digi_particles();
    
    for (size_t j=0;j < parts.size();j++){
      countification = 0;  
      std::vector<bhep::hit*> hits = parts[j]->hits("tracking");
      if ( hits.size() == 0) continue;
      nhit1 = (int)hits.size();
      nhit3 = 0;
      //std::cout << "Voxels before: " << nhit1 << std::endl;
      for (int k=0;k<nhit1;k++){
	x1[k] = hits[k]->x()[0];
	y1[k] = hits[k]->x()[1];
	vdouble dX = hits[k]->vddata("Xpoint");
	vdouble dY = hits[k]->vddata("Ypoint");
	vdouble dZ = hits[k]->vddata("Zpoint");
	vdouble dE = hits[k]->vddata("Epoint");
	//nhit3 += (int)dX.size();
	rawE[k] = 0;
	for (int vdiff=0;vdiff<(int)dX.size();vdiff++){
	  voxXres->Fill( x1[k]-dX[vdiff] );
	  voxYres->Fill( y1[k]-dY[vdiff] );
	  xR[nhit3] = dX[vdiff];
	  yR[nhit3] = dY[vdiff];
	  zR[nhit3] = dZ[vdiff];
	  nhit3++;
	  rawE[k] += dE[vdiff];
	}
	z1[k] = hits[k]->x()[2];
	vE[k] = hits[k]->ddata("TotalEng");
	xVe[0][k] = hits[k]->ddata("XEng");
	xVe[1][k] = hits[k]->ddata("YEng");
      }
      
      hct->execute(hits, meas);
      nhit2 = (int)meas.size();
      for (int l=0;l<nhit2;l++){
	x2[l] = meas[l]->position()[0];
	y2[l] = meas[l]->position()[1];
	z2[l] = meas[l]->position()[2];
	vector<bhep::hit*> vs = meas[l]->get_hits();
	double xx=0, yy=0,QQ=0,Q1;
	nVox[l] = (int)meas[l]->get_nVox();
	nvoxX[l] = meas[l]->get_VoxX();
	cout << meas[l]->get_VoxX() << endl;
	nvoxY[l] = meas[l]->get_VoxY();
	for (size_t ll=0;ll<vs.size();ll++){
	  vdouble XX = vs[ll]->vddata("Xpoint");
	  vdouble YY = vs[ll]->vddata("Ypoint");
	  vdouble EE = vs[ll]->vddata("Epoint");
	  for (size_t lll=0;lll<XX.size();lll++){
	    Q1 = EE[lll];//vs[ll]->ddata( "truEng" );
	    QQ += Q1;
	    xx += XX[lll]*Q1;//vs[ll]->x()[0]*Q1;
	    yy += YY[lll]*Q1;//vs[ll]->x()[1]*Q1;
	  }
	  XX.clear(); YY.clear(); EE.clear();
	}
	//cout << "Weight mean: " << xx/QQ << endl;
	xRes[l] = x2[l] - xx/QQ;
	yRes[l] = y2[l] - yy/QQ;
      }
      // for (int ii=0;ii<nhit2;ii++){
// 	countification += meas[ii]->get_nVox();
// 	std::cout << "MuProp meas["<<ii<<"] = " << meas[ii]->get_mu_prop() << std::endl;
//       }
      
      evt = i;
      tree->Fill();
    }
    
    stc_tools::destroy(meas);
    parts.clear();
    e.clear();

  }

  inDst.close();
  
  f1->Write();
  f1->Close();

  delete hct;

  return 0;
}

/******************************************************
 *                                                    *
 * Class definition of the hit_reconstructor which    *
 * will create the rec_hit's from the deposits        *
 * from the simulation.                               *
 *                                                    *
 * author: Andrew Laing, 11/2009                      *
 *                                                    *
 ******************************************************/

#ifndef _GDML_HIT_CONSTRUCTOR__
#define _GDML_HIT_CONSTRUCTOR__

//#include <recpack/Measurement.h>

#include <bhep/gstore.h>
#include <bhep/hit.h>
#include <bhep/particle.h>

// #include <CLHEP/Random/RanluxEngine.h>
#include <TRandom3.h>

//#include "TFile.h"
#include "TH1F.h"

using namespace bhep;

//#include <mind/rec_hit.h>

#include <map>

class gdml_hit_constructor
{
 public:
  //constructor
  gdml_hit_constructor(const bhep::gstore& store);
  //destructor
  ~gdml_hit_constructor();

  //reconstruction.
  void execute(const std::vector<bhep::hit*>& hits, std::vector<bhep::hit*>& rec_hit, std::vector<TH1F*>& histo_vec);

  //Temporary histograms used for debugging.
  TH1F* rawHitsTH1F;
  TH1F* clusteredHitsTH1F;
  TH1F* digitizedHitsTH1F;
  TH1F* xeTH1F;
  TH1F* xeAttTH1F;
  TH1F* xeSmearTH1F;
  TH1F* yeTH1F;
  TH1F* yeAttTH1F;
  TH1F* yeSmearTH1F;

 private:

  //reset
  void reset();

  //clusteres the hits in xy to get the better resolution expected by the overlaying bars
  void clustering(const std::vector<bhep::hit*>& sortedHits);
  void clustering2(const std::vector<bhep::hit*>& sortedHits);
  void clusteringXY(const std::vector<bhep::hit*> hits, int key);

  void ClusteringHits(const std::vector<bhep::hit*> hits, int key);
  void ClusteringHits2(const std::vector<bhep::hit*> hits, int key);
  std::vector<bhep::hit*> FilteringBadHits(const std::vector<bhep::hit*> hits);

  std::vector<std::vector<bhep::hit*> > NextCluster(std::vector<bhep::hit*> hits,double tolerance, string data);
  //void NextCluster(std::vector<bhep::hit*> hits,double tolerance, string data,std::vector<std::vector<bhep::hit*> > hitsVector);

  bool CorrectHit(const std::vector<bhep::hit*> hits);

  int calculate_vox_no(std::vector<bhep::hit*> hits);

  //Make the rec hits.
  void construct_hits(std::vector<bhep::hit*>& rec_hit);
  //make an individual rec hit.
  bhep::hit* get_vhit(int vox, double z, const std::multimap<int,bhep::hit*>& map1);
  
  //Random Generator for the smear.
  TRandom3 _ranGen;

  //Parameters related to current MIND setup.
  double _detectorLength;
  double _detectorX;
  double _detectorY;
  double _vertexDetdepth;
  //double _vertexDetX;
  //double _vertexDetY;
  //double _passiveLength;
  double _activeLength;
  //double _braceLength;
  //double _gapLength;
  //int _nActive;
  //int OctGeom;
  double _minEng;
  double _attLength;
  
  //EMatrix _cov;

  string _measType;

  //Voxel properties.
  //double _voxXdim;
  //double _voxYdim;
  int _nVoxX;
  //int _nVox;

  //Container wchich will define voxels.
  std::map<double, std::multimap<int, bhep::hit*> > _voxels;
  
};

class forwardSortX{
 public:
  bool operator()(const bhep::hit* p1, const bhep::hit* p2){
    if (p2->x()[0] > p1->x()[0]) return true;
    return false;
  }

};

class forwardSortY{
 public:
  bool operator()(const bhep::hit* p1, const bhep::hit* p2){
    if (p2->x()[1] > p1->x()[1]) return true;
    return false;
  }

};

class forwardSortZ{
 public:
  bool operator()(const bhep::hit* p1, const bhep::hit* p2){
    if (p2->x()[2] > p1->x()[2]) return true;
    return false;
  }

};

class timeSort{
 public:
  bool operator()(bhep::hit* p1, bhep::hit* p2){
    if (p2->ddata( "time" ) > p1->ddata( "time" )) return true;
    return false;
  }

};

#endif

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
  void execute(const std::vector<bhep::hit*>& hits, std::vector<bhep::hit*>& rec_hit);

 private:

  //reset
  //void reset();

  //calculate z position of layers.
  void calculate_layerZ();
  //function which puts hits into the map.
  void parse_to_map(const std::vector<bhep::hit*> hits);
  //find plane of hit.
  double find_plane(bhep::hit& curHit);
  //find vox number of hit.
  int calculate_vox_no(bhep::hit& curHit);
  //Make the rec hits.
  void construct_hits(const std::vector<bhep::hit*>& hits,std::vector<bhep::hit*>& rec_hit);
  //make an individual rec hit.
  bhep::hit* get_vhit(bhep::hit* curr_hit);
  
  //Random Generator for the smear.
  TRandom3 _ranGen;

  //vector of plane z positions.
  std::vector<double> _zLayer;
  //iterator for z planes.
  std::vector<double>::iterator _zIt;

  //Parameters related to current MIND setup.
  double _detectorLength;
  double _detectorX;
  double _detectorY;
  double _vertexDetdepth;
  double _vertexDetX;
  double _vertexDetY;
  double _passiveLength;
  double _activeLength;
  double _braceLength;
  double _gapLength;
  int _nActive;
  int OctGeom;
  double _minEng;
  double _attLength;
  
  //EMatrix _cov;

  string _measType;

  //Voxel properties.
  double _voxXdim;
  double _voxYdim;
  int _nVoxX;
  int _nVox;

  //Container wchich will define voxels.
  //std::map<double, std::multimap<int, bhep::hit*> > _voxels;
  
};

class forwardSort{
 public:
  bool operator()(const bhep::hit* p1, const bhep::hit* p2){
    if (p2->x()[2] > p1->x()[2]) return true;
    return false;
  }

};

#endif

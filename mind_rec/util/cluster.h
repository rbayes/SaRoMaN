/******************************************************
 *                                                    *
 * Class definition of cluster for hits formed        *
 * from many reconstructed hits in MIND               *
 *                                                    *
 * author: Andrew Laing, 11/2009                      *
 *                                                    *
 ******************************************************/

#ifndef _CLUSTER__
#define _CLUSTER__

#include <bhep/hit.h>
#include <bhep/gstore.h>

#include <recpack/Measurement.h>

class cluster: public Measurement
{
 public:
  //Constructor.
  cluster();
  //destructor.
  ~cluster();

  //add a contibuting hit.
  void add_hit(bhep::hit* dep);

  //Getters.
  double get_eng(){ return _eng; }
  int get_nhits(){ return _nhit; }
  size_t get_nVox(){ return _nVox; }
  int get_VoxX(){ return _nVoxV[0]; }
  int get_VoxY(){ return _nVoxV[1]; }
  double get_mu_prop(){ return _muProp / (double)_nhit; }
  std::vector<bhep::hit*> get_hits(){ return _voxes; }
  double get_time(){return _time;}

  bhep::vstring get_mother_particle(){return _mother_particle;}

  //Setters.
  void set_VoxX(int n){ _nVoxV[0] = n; }
  void set_VoxY(int n){ _nVoxV[1] = n; }
  
 private:

  std::vector<bhep::hit*> _voxes;

  size_t _nVox;

  int _nhit;

  int _nVoxV[2];

  double _eng;

  double _muProp;

  double _time;

  bhep::vstring _mother_particle;

};

#endif

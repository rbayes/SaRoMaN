/******************************************************
 *                                                    *
 * Class definition of rec_hit for hits reconstructed *
 * from many particle an shower hits in MIND          *
 *                                                    *
 * author: Andrew Laing, 11/2009                      *
 *                                                    *
 ******************************************************/

#ifndef _rec_hit__
#define _rec_hit__

#include <recpack/Measurement.h>

#include <bhep/hit.h>

class rec_hit: public Measurement
{
 public:
  //Constructor
  rec_hit();
  //Destructor
  ~rec_hit();

  //Add a point.
  void add_point(bhep::hit* point);

  //Get total hit energy.
  double get_hit_eng() { return _totEng; }
  //Get no. contributers.
  unsigned int get_npoints() { return _npoints; }
  //Get muon proportion of contributers.
  double get_mu_prop() { return _muProp/(double)_npoints; }
  //Get vector of contributers.
  std::vector<bhep::hit*> get_points() { return _points; }
  //set the vox number.
  void set_vox_no(int num){ _voxNo = num; }
  //get the vox number.
  int get_vox_no(){ return _voxNo; }

 private:

  std::vector<bhep::hit*> _points;//vector of contributing deposits.

  unsigned int _npoints;//number of contributions

  double _totEng;

  double _muProp;//measure of 'muon-ness'

  int _voxNo;

};

#endif

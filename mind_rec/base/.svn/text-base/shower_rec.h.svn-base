/* -*- mode: c++ -*- */
#ifndef _shower_rec__
#define _shower_rec__

#include <recpack/Measurement.h>

#include <TGraphErrors.h>
#include <TF1.h>

#include <mind/shower.h>

#define NMAXLAYERS 150

class shower_rec{

 public:

  shower_rec();

  ~shower_rec(){};

  //Main control functions.
  bool initialize( double tol, double resolution );
  bool execute( shower& had_shower );
  bool finalize();
  //

 private:

  void order_hits( measurement_vector& hads );
  bool reconstruct( EVector vert );
  void reset();

  double _tolerance;
  double _res;
  EVector _dirVec;
  measurement_vector _meas_in_plane[NMAXLAYERS];
  std::vector<double> _zpos;
  size_t _planes;

  //protected:

};

class forwardSorter{
public:
  bool operator()(const Measurement* p1, const Measurement* p2){
    if (p2->position()[2] > p1->position()[2]) return true;
    return false;
  }
  
};

#endif

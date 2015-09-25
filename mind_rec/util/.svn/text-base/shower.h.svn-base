/* -*- mode: c++ -*- */
#ifndef _shower__
#define _shower__

#include <recpack/Measurement.h>

class shower{

 public:

  shower( const measurement_vector& hads, EVector& vert_loc );

  ~shower(){};

  //functions to set values
  //direction unit vector.
  void set_unit_direction( EVector dir )
    { _direction = dir; }
  //Four momentum.
  void set_four_mom( EVector mom )
    { _fourMom = mom; }

  //Getters.
  //direction.
  EVector& get_direction()
    { return _direction; }
  //Four momentum.
  EVector& get_fourMom()
    { return _fourMom; }
  //Vertex.
  EVector& get_vertex()
    { return _vertex; }
  //3 momentum.
  EVector get_threeMom()
    { return Recpack::box( _fourMom, 0, 3 ); }
  //Energy.
  double get_energy()
    {return _fourMom[3]; }
  //All hits.
  measurement_vector& get_had_hits()
    { return _hits; }

 private:

  EVector _direction;
  EVector _vertex;
  EVector _fourMom;

  measurement_vector _hits;

};

#endif

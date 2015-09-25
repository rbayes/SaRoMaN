#include <shower.h>

shower::shower( const measurement_vector& hads, EVector& vert_loc )
{

  _vertex = vert_loc;
  _direction = EVector(3,0);
  _fourMom = EVector(4,0);

  measurement_vector::const_iterator It1;
  for (It1 = hads.begin();It1 != hads.end();It1++)
    _hits.push_back( (*It1) );

}

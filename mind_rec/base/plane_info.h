#ifndef _plane_info___
#define _plane_info___


#include <mind/cluster.h>
#include <mind/MINDsetup.h>///without this bhep::unit not working ?
#include <bhep/gstore.h>


class plane_info{
  
 public:
  
  
  //constructor
  // plane_info(int no, double z, const bhep::gstore& pstore);
  plane_info();
  plane_info(int no, double z, const bhep::gstore& pstore);
  plane_info(int no, double z, double r, const bhep::gstore& pstore);
  
  ///destructor
  virtual  ~plane_info();
    
  ///get functions
  std::vector<cluster*>& GetHits() {return _hits;}
  int GetNHits() {return _hits.size();}
  double GetZ()  {return _zPos ;}
  double GetR()  {return _rPos ;}
  double GetEng() { return _engPlane;}
  int GetPlaneNo() { return _planeNo;}


  ///add hit to the plane
  void AddHit(cluster* hit);

  ///energy deposition
  double correctEdep(double edep, cluster* hit);

  void reset();
  
  //Return index of outlying hit
  int mostSeparatedHitIndex();

 protected:

  
  ///data members

  std::vector<cluster*> _hits;
  double _zPos;
  double _engPlane;
  int _planeNo;
  double _rPos;
  bhep::gstore _store;

 
};
#endif

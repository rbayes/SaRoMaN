/******************************************************
 *                                                    *
 * Class definition of the hit_clusterer which        *
 * will create the clusters from the rec_hits         *
 * from the digitization.                             *
 *                                                    *
 * author: Andrew Laing, 11/2009                      *
 *                                                    *
 ******************************************************/

#ifndef _HIT_CLUSTERER__
#define _HIT_CLUSTERER__

#include <mind/cluster.h>

#include <bhep/hit.h>

#include <CLHEP/Random/RanluxEngine.h>

class hit_clusterer
{
 public:
  //constructor.
  hit_clusterer(const bhep::gstore& store);
  //destructor.
  ~hit_clusterer();

  void execute(const std::vector<bhep::hit*>& deps,
	       std::vector<cluster*>& clusts);

 private:

  //make clusters in a particular plane;
  void clusters_in_plane(double zpos, std::map<int,bhep::hit*>& deps,
			 std::vector<cluster*>& clusts);

  void form_clusters(double zpos, std::map<int,bhep::hit*>& deps,
			 std::vector<cluster*>& clusts);

  void calculate_clust_pos(const std::vector<bhep::hit*>& hits, EVector& vec);

  //make cluster from one hit.
  cluster* make_cluster(double zpos, bhep::hit* dep);
  //make cluster from many hits.
  cluster* make_cluster(const EVector& vec,
			const std::vector<bhep::hit*>& deps);

  //get the vox with maximum energy.
  int get_max_vox(const std::map<int,double>& voxes);

  //Random Generator for the smear.
  //RanluxEngine _ranGen;
  //Voxels in X and Y.
  int _nVoxV[2];

  //sigma for position smear.
  double _sigMa;
  double _sigMaZ;
  double _minEng;
  double _res[3];

  //number of voxels per x edge.
  int _vInX;

  //covarience and measurement type.
  EMatrix _cov;
  string _measType;

};

#endif

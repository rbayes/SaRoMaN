#ifndef _super_fit___
#define _superfit___

#include <recpack/RecPackManager.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/LsqFitter.h>
#include <recpack/ParticleState.h>

#include <mind/MINDsetup.h>
#include <mind/Utilities.h>
#include <mind/MINDfitman.h>
#include <mind/cluster.h>
//#include <mind/event_classif.h>
#include <mind/hit_clusterer.h>

#include <bhep/event.h>
#include <bhep/gstore.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>

using namespace Recpack;

class super_fit{

 public:

  super_fit();
  
  virtual ~super_fit();

  MINDsetup _supergeom;

  // Calculates momenta using a simple range calculation.
  double RangeMomentum(double length,double nodeZ);
  // Calculates momenta from the curvature of the track, should be used for momenta greater
  // than can be calculated from RangeMomentum
  double MomentumFromCurvature(const Trajectory& traj, int startPoint,double minP);
  // Calculate the charge by using a quadratic fit.
  double CalculateCharge(const Trajectory& startTrack);

};

#endif

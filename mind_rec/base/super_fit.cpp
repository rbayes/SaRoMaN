#include <super_fit.h>
#include <TMath.h>

using namespace bhep;///

//**********************************************************************
super_fit::super_fit() {
  //**********************************************************************
}

//***********************************************************************
super_fit::~super_fit() {
  //***********************************************************************
}


/***************************************************************************************/
double super_fit::RangeMomentum(double length,double nodeZ){
  /***************************************************************************************/
  // Get momentum depending on length. Did it go through each submodule? Get specific wFe
  // for each submodule.

  std::map<dict::Key,vector<double> > moduleDataMap = _supergeom.getModuleDataMap();
  double p = 0;

  for (std::map<dict::Key,vector<double> >::iterator it=moduleDataMap.begin();
       it!=moduleDataMap.end(); ++it)
    {
      double module_pos = it->second[0];
      double module_half_size = it->second[1];
      double wFe = it->second[2];

      //Sanity checking not outside range forward 
      if(!((module_pos-module_half_size)<(nodeZ+length)))
	{ 
	  continue;
	}

      // Through the whole module
      else if((nodeZ + length) > (module_pos + module_half_size))
	{
	   p += 2*module_half_size * wFe;
	}
      // Stop in the module
      else if((nodeZ + length) < (module_pos + module_half_size))
	{
	  p+=((length + nodeZ) - (module_pos - module_half_size))*wFe;
	}	 
    }

  return p;

}

/***************************************************************************************/
double super_fit::MomentumFromCurvature(const Trajectory& traj, int startPoint,double minP){
  /***************************************************************************************/

  // Unsure which hit is first.

  bool startLow = false;

  if(traj.nodes()[startPoint]->measurement().vector()[2] < traj.nodes()[traj.size()-1]->measurement().vector()[2])
    {
      startLow = true;
    }

  // Find straight line between 2 and 3 point.
  // y = k*z+m

  int start = startPoint;
  int point = traj.size()-3;

  double scalar = dot(_supergeom.getRawBField(traj.nodes()[start]->measurement().vector()),
		     _supergeom.getRawBField(traj.nodes()[point+2]->measurement().vector()));

  while(scalar < 0)
    {
      if(startLow)
	{
	  point--;
	}
      else
	{
	  start++;
	}
      scalar = dot(_supergeom.getRawBField(traj.nodes()[start]->measurement().vector()),
		     _supergeom.getRawBField(traj.nodes()[point+2]->measurement().vector()));
    }


  double k1 = 0;

  while(k1 == 0 && start <= point)
    {
      k1 = ((traj.nodes()[start]->measurement().vector()[1]-traj.nodes()[start+1]->measurement().vector()[1])/
	    (traj.nodes()[start]->measurement().vector()[2]-traj.nodes()[start+1]->measurement().vector()[2]) +
	    (traj.nodes()[start+1]->measurement().vector()[1]-traj.nodes()[start+2]->measurement().vector()[1])/
	    (traj.nodes()[start+1]->measurement().vector()[2]-traj.nodes()[start+2]->measurement().vector()[2]))/2;
      
      start++;
    }

  double k2 = 0;

  while(k2 == 0 & 0 <= point)
    {
      k2 = ((traj.nodes()[point]->measurement().vector()[1]-traj.nodes()[point+1]->measurement().vector()[1])/
	    (traj.nodes()[point]->measurement().vector()[2]-traj.nodes()[point+1]->measurement().vector()[2]) +
	    (traj.nodes()[point+1]->measurement().vector()[1]-traj.nodes()[point+2]->measurement().vector()[1])/
	    (traj.nodes()[point+1]->measurement().vector()[2]-traj.nodes()[point+2]->measurement().vector()[2]))/2;

      point--;
    }

  double m1 = traj.nodes()[1]->measurement().vector()[1] - k1* (traj.nodes()[1]->measurement().vector()[2]);  
  double m2 = traj.nodes()[point+1]->measurement().vector()[1] - k2* (traj.nodes()[point+1]->measurement().vector()[2]);
    
  // Find their orthogonal lines.

  double k3 = -1/k1;
  double k4 = -1/k2;

  double m3 = traj.nodes()[1]->measurement().vector()[1] - k3* (traj.nodes()[1]->measurement().vector()[2]); 
  double m4 = traj.nodes()[point+1]->measurement().vector()[1] - k4* (traj.nodes()[point+1]->measurement().vector()[2]); 

  // Find intersection

  double zc = (m3-m4)/(k4-k3);

  double yc = k3 *zc + m3;

  double r = sqrt((zc-traj.nodes()[1]->measurement().vector()[2])*(zc-traj.nodes()[1]->measurement().vector()[2])
		  + (yc-traj.nodes()[1]->measurement().vector()[1])*(yc-traj.nodes()[1]->measurement().vector()[1]));

  double r2 = sqrt((zc-traj.nodes()[point+1]->measurement().vector()[2])*(zc-traj.nodes()[point+1]->measurement().vector()[2])
		  + (yc-traj.nodes()[point+1]->measurement().vector()[1])*(yc-traj.nodes()[point+1]->measurement().vector()[1]));

  //EVector B = _supergeom.getBField(pos);

  double field = _supergeom.getRawBField(traj.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr();

  double correction = 3.5;

  double finalP = 300*fabs(field) * correction *(r+r2)/2;

  cout<<"Momentum highPt final super"<<finalP<<endl;

  if(finalP<minP)
    {
      // Try to recursivly find a better value.
      if((start +1)< point)
	{
	  finalP =  MomentumFromCurvature(traj,start+1,minP);
	}
      else
	{
	  finalP = minP;
	}
    }
  // Need a better estimate of a maximum value for the detector.
  if(finalP > 14000)
    {
      // Try to recursivly find a better value.
      if((start +1)< point)
	{
	  finalP =  MomentumFromCurvature(traj,start+1,minP);
	}
      else
	{
	  finalP = 14000;
	}
    }

  return finalP;
}

//***********************************************************************
double super_fit::CalculateCharge(const Trajectory& track) {
  //***********************************************************************

  int fitcatcher;
  int nMeas = track.size();

  // Run this with different track starts! Will make really bad guesses when we change field.
  // Either from start or from the end.

  //Try both start from back, only 4 hits when pt to large?

  //if (nMeas > 4) nMeas = 4;

  // If we pass through the detector, remove the last 2? Plane hits.

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;

  double firstNodeZ = track.nodes()[nMeas-1]->measurement().position()[2];
  
  pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
  pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
  pos[2] = z[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[2];

  EVector prevB(3,0);

  for (int iMeas = nMeas-2;iMeas >= 0;iMeas--){
    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];

    // get the b-field from the previous step
    EVector B = _supergeom.getBField(pos);

    if(dot(prevB,B) < 0)
      {
	break;
      }
    prevB=B;

    pos[0] = x[iMeas];  pos[1] = y[iMeas];   pos[2] = z[iMeas];

    if(crossprod(Z,B).norm() ==0 )
      {
	u[iMeas] = 0;
      }
    else
      {
	u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
      }
  }

  TGraph *gr = new TGraph((const int)minindex, z, u);
  
  TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun->SetParameters(0.,0.001,0.001);

  fitcatcher = gr->Fit("parfit2", "QN");

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -

  double qtilde = 0.3*pow(1+pow(fun->GetParameter(1),2),3./2.)/
    (2*fun->GetParameter(2));

  int meansign = (int)(qtilde/fabs(qtilde));
  

  // Better to be wrong than to return 0.
  if(meansign == 0)
    {
      meansign = 1;
    }

  return meansign;
}

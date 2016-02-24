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

      //std::cout<<"Fitter "<<module_pos<<" "<<module_half_size<<" "
      //     <<wFe<<" "<<nodeZ<<" "<<length<<std::endl;

      std::cout<<"Fitter "<<nodeZ<<" "<<length<<" "<<wFe<<" "<<module_pos<<std::endl;

      //Sanity checking not outside range forward 
      if(!((module_pos-module_half_size)<(nodeZ+length)))
	{ 
	  continue;
	  //std::cout<<"Fitter not in range"<<std::endl;
	  
	}

      // Through the whole module
      else if((nodeZ + length) > (module_pos + module_half_size))
	{
	   p += 2*module_half_size * wFe;
	   std::cout<<"p+="<<2*module_half_size * wFe<<std::endl;
	   std::cout<<"Fitter through"<<std::endl;
	}
      // Stop in the module
      else if((nodeZ + length) < (module_pos + module_half_size))
	{
	  p+=((length + nodeZ) - (module_pos - module_half_size))*wFe;
	  //p+=(module_half_size + module_pos -(length+nodeZ))*wFe;
	  std::cout<<"p+="<<((length + nodeZ) - (module_pos - module_half_size))*wFe<<std::endl;
	  std::cout<<"Fitter stop"<<std::endl;
	}	 
    }

  //std::cout<<"In Map p "<<p<<std::endl;

  return p;

}

/***************************************************************************************/
double super_fit::MomentumFromCurvature2(const Trajectory& traj){
  /***************************************************************************************/

  // Calculate the radius of curvature by using the available hits from the same magnetic 
  // region (same sign) and same submodule (same magnitude)

  std::map<dict::Key,vector<double> > moduleDataMap = _supergeom.getModuleDataMap();
  vector<vector<EVector> > moduleSplitHits;

  //Group hits in each submodule
  for (std::map<dict::Key,vector<double> >::iterator it=moduleDataMap.begin();
       it!=moduleDataMap.end(); ++it)
    {
      vector<EVector> moduleHits;

      double module_pos = it->second[0];
      double module_half_size = it->second[1];

      EVector lastB(3,0);

      for(int cnt = 0; cnt < traj.size(); cnt++)
	{
	  EVector pos(3,0);
	  pos[0] = traj.nodes()[cnt]->measurement().vector()[0];
	  pos[1] = traj.nodes()[cnt]->measurement().vector()[1];
	  pos[2] = traj.nodes()[cnt]->measurement().vector()[2];
	  //cout<<"In MomentumFromDeflection: "<<pos[2]<<endl; 
	  //cout<<module_pos<<" "<<module_half_size<<endl;

	  // Not in the module?
	  if(!((module_pos - module_half_size)<pos[2] && pos[2]<(module_pos + module_half_size)))
	    {
	      // Jump to next value in loop.
	      //cout<<"jump"<<endl;
	      continue;
	    }

	  EVector B = _supergeom.getBField(pos);

	  if(dot(lastB,B)< 0)
	    {
	      // B changed sign, break first loop.
	      //cout<<"bchange"<<endl;
	      break;
	    }
	  
	  lastB = B;
	  //cout<<"fill"<<endl;
	  moduleHits.push_back(pos);

	}//End for traj

      //cout<<"fill"<<endl;
      moduleSplitHits.push_back(moduleHits);
	 
    }//End for moduleDataMap

  double finalP = 0;
  double applicablemodules =0;

  for(int cnt = 0; cnt < moduleSplitHits.size(); cnt++)
    {
      double moduleP = 0;
      double applicable =0;

      //if((moduleSplitHits[cnt].size()-1) < 3)
      if((moduleSplitHits[cnt].size()) < 5)
	{
	  continue;
	}

      //double y[(const int) moduleSplitHits[cnt].size()], z[(const int) moduleSplitHits[cnt].size()];

      applicablemodules++;

      const int nMeas = moduleSplitHits[cnt].size()-1;

      double x[(const int)nMeas], y[(const int)nMeas], 
	z[(const int)nMeas];
      
      for(int icnt = 0; icnt < moduleSplitHits[cnt].size()-1; icnt++)
	{
	  y[icnt] = moduleSplitHits[cnt][icnt+0][1];
	  z[icnt] = moduleSplitHits[cnt][icnt+0][2];
	}
      
      TCanvas *c1 = new TCanvas();

      TGraph *gr = new TGraph((const int)nMeas, z, y);
      
      TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
      fun->SetParameters(0.,0.001,0.001);
      
      int fitcatcher = gr->Fit("parfit2", "QN");

      gr->Draw();

      c1->SaveAs("here.pdf");

      double qtilde = 0.3*pow(1+pow(fun->GetParameter(1),2),3./2.)/
	(2*fun->GetParameter(2));      


      //double qtilde = 0.3*pow(1+pow(fun->GetParameter(1)+2*fun->GetParameter(2)*z[0],2),3./2.)/
      //(2*fun->GetParameter(2));

      cout<<"New P "<<qtilde<<endl;
      moduleP += qtilde; // instead of the for loop. 
      applicable++;

      /*
      for(int icnt = 0; icnt < (moduleSplitHits[cnt].size() -2); icnt++)
	{
	  cout<<"module: "<<cnt<<endl;
	  cout<<"hit: "<<icnt<<endl;

	 
	  //cout<<"Y2: "<<moduleSplitHits[cnt][icnt+2][1]<<endl;
	  //cout<<"Z2: "<<moduleSplitHits[cnt][icnt+2][2]<<endl;

	  EVector pos2(2,0); //(x,y,z)
	  //pos2[0] = moduleSplitHits[cnt][icnt][0];
	  pos2[0] = moduleSplitHits[cnt][icnt+0][1];
	  pos2[1] = moduleSplitHits[cnt][icnt+0][2];
	  
	  EVector pos1(2,0);
	  //pos1[0] = moduleSplitHits[cnt][icnt+1][0];
	  pos1[0] = moduleSplitHits[cnt][icnt+1][1];
	  pos1[1] = moduleSplitHits[cnt][icnt+1][2];

	  EVector pos0(2,0);
	  //pos0[0] = moduleSplitHits[cnt][icnt+2][0];
	  pos0[0] = moduleSplitHits[cnt][icnt+2][1];
	  pos0[1] = moduleSplitHits[cnt][icnt+2][2];


	  // See explination in documentation or on t.ex.: 
	  // http://www.regentsprep.org/regents/math/geometry/gcg6/RCir.htm

	  double y3 = moduleSplitHits[cnt][icnt+0][1];
	  double y2 = moduleSplitHits[cnt][icnt+1][1];
	  double y1 = moduleSplitHits[cnt][icnt+2][1];

	  double z3 = moduleSplitHits[cnt][icnt+0][2];
	  double z2 = moduleSplitHits[cnt][icnt+1][2];
	  double z1 = moduleSplitHits[cnt][icnt+2][2];

	  if(z3<z1)
	    {
	      z3 = moduleSplitHits[cnt][icnt+2][2];
	      y3 = moduleSplitHits[cnt][icnt+2][1];

	      z1 = moduleSplitHits[cnt][icnt+0][2];
	      y1 = moduleSplitHits[cnt][icnt+0][1];
	    }

	  cout<<"Z1: "<<z1<<endl;
	  
	  cout<<"Z2: "<<z2<<endl;

	  cout<<"Z3: "<<z3<<endl;


	  double mr = (y2-y1)/(z2-z1);
	  double mt = (y3-y2)/(z3-z2);

	  double zc = (mr*mt*(y3-y1)+mr*(z2+z3)-mt*(z1+z2))/(2*(mr-mt)); 

	  double yc = -1/mr * (zc - (z1+z2)/2) + (y1+y2)/2;

	  double radius = 1000 * sqrt((y1-yc)*(y1-yc) + (z1-zc)*(z1-zc));

	  cout<<"Radius: "<<radius<<endl;

	  double B =  fabs(_supergeom.getBField(moduleSplitHits[cnt][icnt+2])[0]); //Should do a cross to get B field.
      
	  double p = 0.3* B * radius;

	  cout<<"Momentum highPt super "<<p<<endl;
	  
	  p = fabs(p);
	  

	  //p==p check for nan.
	  // did we calculate a resonable p?
	  if( p< 8000 &&!(p != p))
	    {
	      moduleP += p; 
	      applicable++;
	    }
	  else
	    {
	      p = 0;
	    }
	}
      */


      //moduleP = moduleP/(moduleSplitHits[cnt].size() -2);
      if(moduleP != 0)
	{
	  moduleP = moduleP/applicable;
	}
      else
	{
	  moduleP = 0;
	}   
      cout<<"Momentum highPt module super"<<moduleP<<endl;

      finalP +=moduleP;
}

  if(finalP != 0)
    {
      finalP = finalP/applicablemodules;
    }
  else
    {
      finalP = 0;
    }
  
  cout<<"Momentum highPt final super"<<finalP<<endl;


  return finalP;
}

/***************************************************************************************/
double super_fit::MomentumFromCurvature(const Trajectory& traj, int startPoint){
  /***************************************************************************************/

  // Unsure which hit is first.

  bool startLow = false;

  if(traj.nodes()[0]->measurement().vector()[2] < traj.nodes()[traj.size()-1]->measurement().vector()[2])
    {
      startLow = true;
    }

  // Find straight line between 2 and 3 point.
  // y = k*z+m

  int start = startPoint;
  int point = traj.size()-3;
  /*
  int final = traj.size()-1;

  //new test using the geometry
  if(startLow)
    {
      k1 = (traj.nodes()[3]->measurement().vector()[1]-traj.nodes()[2]->measurement().vector()[1])/
	(traj.nodes()[3]->measurement().vector()[2]-traj.nodes()[2]->measurement().vector()[2]);

      k2 = (traj.nodes()[final]->measurement().vector()[1]-traj.nodes()[final-1]->measurement().vector()[1])/
	    (traj.nodes()[final]->measurement().vector()[2]-traj.nodes()[final-1]->measurement().vector()[2]);
    }
  else
    {
      k1 = (traj.nodes()[final-3]->measurement().vector()[1]-traj.nodes()[final-2]->measurement().vector()[1])/
	(traj.nodes()[final-3]->measurement().vector()[2]-traj.nodes()[final-2]->measurement().vector()[2]);

      k2 = (traj.nodes()[0]->measurement().vector()[1]-traj.nodes()[1]->measurement().vector()[1])/
	    (traj.nodes()[0]->measurement().vector()[2]-traj.nodes()[1]->measurement().vector()[2]);
    }
  */

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

  cout<<"line 1: "<<k1<<"z+"<<m1<<endl;
  cout<<"line 2: "<<k2<<"z+"<<m2<<endl;
  cout<<"line 3: "<<k3<<"z+"<<m3<<endl;
  cout<<"line 4: "<<k4<<"z+"<<m4<<endl;

  double zc = (m3-m4)/(k4-k3);

  double yc = k3 *zc + m3;

  cout<<"yc: "<<yc<<" zc: "<<zc<<endl;


  double r = sqrt((zc-traj.nodes()[1]->measurement().vector()[2])*(zc-traj.nodes()[1]->measurement().vector()[2])
		  + (yc-traj.nodes()[1]->measurement().vector()[1])*(yc-traj.nodes()[1]->measurement().vector()[1]));

  double r2 = sqrt((zc-traj.nodes()[point+1]->measurement().vector()[2])*(zc-traj.nodes()[point+1]->measurement().vector()[2])
		  + (yc-traj.nodes()[point+1]->measurement().vector()[1])*(yc-traj.nodes()[point+1]->measurement().vector()[1]));

  cout<<"Final r: "<<r<<endl;
  cout<<"Final r2: "<<r2<<endl;

  // Need better value for b;

  // 1.25 correction factor.


  //EVector B = _supergeom.getBField(pos);

  double field = _supergeom.getRawBField(traj.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr();

  double correction = 3.5;

  //cout<<"Scale "<<_supergeom.getBScaleAvr()<<endl;
  //cout<<"field "<< 3.5*_supergeom.getRawBField(traj.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr()<<endl;
  //cout<<"current "<<0.000957447*1.25<<endl;
  //  double finalP = 300*0.000957447*1.25*(r+r2)/2;


  double finalP = 300*fabs(field) * correction *(r+r2)/2;

  cout<<"Momentum highPt final super"<<finalP<<endl;

  if(finalP<1400)
    {
      // Try to recursivly find a better value.
      if((start +1)< point)
	{
	  finalP =  MomentumFromCurvature(traj,start+1);
	}
      else
	{
	  finalP = 1400;
	}
    }

  if(finalP > 8000)
    {
      // Try to recursivly find a better value.
      if((start +1)< point)
	{
	  finalP =  MomentumFromCurvature(traj,start+1);
	}
      else
	{
	  finalP = 8000;
	}
    }
  
  // if(finalP<1400 || finalP> 8000 || finalP != finalP)
  //  {
  //cout<<"Flag: "<<finalP<<endl;
  //finalP = 0;
  //  }

  return finalP;
}


//***********************************************************************
double super_fit::CalculateCharge(const Trajectory& track) {
  //***********************************************************************

  int fitcatcher;
  int nMeas = track.size();

  //if (nMeas > 4) nMeas = 4;

  // If we pass through the detector, remove the last 2? Plane hits.

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];/// r[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;
  //double pathlength=0;/// tminR=9999.9;

  double firstNodeZ = track.nodes()[nMeas-1]->measurement().position()[2];
  
  pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
  pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
  pos[2] = z[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[2];
  //pos[2] = 0.0;

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
    cout<<"mag B: "<<B[0]<<" "<<B[1]<<" "<<B[2]<<endl;
    if(crossprod(Z,B).norm() ==0 )
      {
	u[iMeas] = 0;
      }
    else
      {
	u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
      }
    cout<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
    cout<<"values: z: "<<z[iMeas]<<" u: "<<u[iMeas]<<endl;
  }

  TGraph *gr = new TGraph((const int)minindex, z, u);
  
  TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun->SetParameters(0.,0.001,0.001);

  fitcatcher = gr->Fit("parfit2", "QN");

  //cout<<"Parameters: "<<fun->GetParameter(1)<<" "<<fun->GetParameter(2)<<endl;

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -

  double qtilde = 0.3*pow(1+pow(fun->GetParameter(1),2),3./2.)/
    (2*fun->GetParameter(2));

  //double qtilde = 0.3*pow(1+pow(fun->GetParameter(1)+2*fun->GetParameter(2)*z[0],2),3./2.)/
  //(2*fun->GetParameter(2));

  cout<<"qtilde: "<<qtilde<<endl;

  int meansign = (int)(qtilde/fabs(qtilde));
  //delete fun;
  //return qtilde;
  //int meansign = sumdq/fabs(sumdq);

  // Handle low qtilde.
  if(fabs(meansign) > 1)
    {
      meansign = 0;
    }
  
  return meansign;
}

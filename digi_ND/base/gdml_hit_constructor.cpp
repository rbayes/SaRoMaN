#include <gdml_hit_constructor.h>

// #include <CLHEP/Random/RandGauss.h>
#include <TRandom3.h>
#include <algorithm>

gdml_hit_constructor::gdml_hit_constructor(const bhep::gstore& store)
{
  //Store detector geometry and calculate scint layer positions.

  _detectorLength = store.fetch_dstore("MIND_z") * m;
  _detectorX = store.fetch_dstore("MIND_x") * m;
  _detectorY = store.fetch_dstore("MIND_y") * m;
  _vertexDetdepth = store.fetch_dstore("vertex_z") * m;
  //_vertexDetX = store.fetch_dstore("vertex_x") * m;
  //_vertexDetY = store.fetch_dstore("vertex_y") * m;
  //_passiveLength = store.fetch_dstore("passive_thickness") * cm;
  _activeLength = store.fetch_dstore("active_thickness") * cm;
  cout<<"_activeLength"<<_activeLength<<endl;
  //_braceLength = 0.0;
  //if (store.find_dstore("bracing_thickness"))
  //  _braceLength = store.fetch_dstore("bracing_thickness") * cm;
  //_gapLength = store.fetch_dstore("air_gap") * cm;
  //_nActive = store.fetch_istore("active_layers");
  //OctGeom = 0;
  //if (store.find_istore("isOctagonal"))
  //  OctGeom = store.fetch_istore("isOctagonal");
  //_voxXdim = store.fetch_dstore("rec_boxX") * cm;
  //_voxYdim = store.fetch_dstore("rec_boxY") * cm;
  //_nVoxX = (int)( _detectorX / _voxXdim );
  //_nVox = _nVoxX * (int)( _detectorY / _voxYdim );
  _nVoxX =store.fetch_istore("nVoxX");

  long seed = (long)store.fetch_dstore("Gen_seed");
  _ranGen = TRandom3( seed );

  _minEng = store.fetch_dstore("min_eng") * MeV;

  _attLength = store.fetch_dstore("WLSatten");  

}

gdml_hit_constructor::~gdml_hit_constructor()
{

}

void gdml_hit_constructor::reset()
{
  //Clear out map.correct??
  _voxels.clear();

}

void gdml_hit_constructor::execute(const std::vector<bhep::hit*>& hits,
				   std::vector<bhep::hit*>& rec_hit, std::vector<TH1F*>& histo_vec)
{
  /*
    Main executible function for the hit_constructor.
  */

  //Set the histogram plotters
  rawHitsTH1F = histo_vec[0];
  clusteredHitsTH1F = histo_vec[1];
  digitizedHitsTH1F = histo_vec[2];
  xeTH1F = histo_vec[3];
  xeAttTH1F = histo_vec[4];
  xeSmearTH1F = histo_vec[5];
  yeTH1F = histo_vec[6];
  yeAttTH1F = histo_vec[7];
  yeSmearTH1F = histo_vec[8];

  //First clear out map.
  reset();
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;
  cout<<"hits size: "<<sortedHits.size()<<endl;

  //sort( sortedHits.begin(), sortedHits.end(), forwardSortX() );
  //sort( sortedHits.begin(), sortedHits.end(), forwardSortY() );
  sort( sortedHits.begin(), sortedHits.end(), forwardSortZ() );

  cout<<"first z: "<<sortedHits[0]->x()[2]<<endl;
  cout<<"last z: "<<sortedHits[sortedHits.size()-1]->x()[2]<<endl;



  clustering(sortedHits);
  //clustering2(sortedHits);

  //Make rec_hits from vox.
  construct_hits( rec_hit );
}

void gdml_hit_constructor::clustering(const std::vector<bhep::hit*>& zSortedHits)
{
  /*
    Cluster the real hits (bar positions from hits) to produce hit positions.
    Also utilize the bar overlap to be able to give an even better position.
    Main jobs is done by calling clusteringXY.
  */  
  std::vector<bhep::hit*>::const_iterator hitIt;
  std::vector<bhep::hit*> moduleHits;
  std::vector<std::vector<bhep::hit*> > moduleHitsVector;
  //std::vector<std::vector<double> > clustered_hits;

  // Fill vectors with hits in the same module (xyxy),sorted by barPosZ.
  for (hitIt = zSortedHits.begin();hitIt != zSortedHits.end();hitIt++)
    {
      double currZ = (*hitIt)->ddata( "barPosZ" );
      double nextZ;
      rawHitsTH1F->Fill((*hitIt)->x()[2]);
      
      if(hitIt + 1 != zSortedHits.end()){ nextZ = (*(hitIt + 1))->ddata( "barPosZ" );}
      else {nextZ = currZ + 3./4. * _activeLength;}

      if(fabs(currZ-nextZ) < 3./4. * _activeLength)
	{
	  moduleHits.push_back((*hitIt));
	}
      else//Next is to far away
	{
	  moduleHits.push_back((*hitIt));

	  moduleHitsVector.push_back(moduleHits);
	  moduleHits.clear();	  
	}    
    }

  //moduleHitsVector per z


  // Retretrive the different hits (if different hits in x or y)

  // Do the actually clustering
  for(int counter = 0; counter < moduleHitsVector.size(); counter++)
    {
      
      cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
	  <<"\t"<<"IsYBar"<<"\t"<<"barNum"
	  <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
	  <<"\t"<<"time"<<"\t"<<"EnergyDep"
	  <<"\t"<<"mother name"<<endl;
      for(int inCounter = 0; inCounter < moduleHitsVector[counter].size(); inCounter++)
	{
	  cout<<moduleHitsVector[counter][inCounter]->x()[0]<<"\t"
	      <<moduleHitsVector[counter][inCounter]->x()[1]<<"\t"
	      <<moduleHitsVector[counter][inCounter]->x()[2]<<"\t"
	      <<moduleHitsVector[counter][inCounter]->idata( "IsYBar" )<<"\t"
	      <<moduleHitsVector[counter][inCounter]->idata( "barNumber" )<<"\t"
	      <<moduleHitsVector[counter][inCounter]->ddata( "barPosZ" )<<"\t"
	      <<moduleHitsVector[counter][inCounter]->ddata( "barPosT" )<<"\t"
	      <<moduleHitsVector[counter][inCounter]->ddata( "time" )<<"\t"
	      <<moduleHitsVector[counter][inCounter]->ddata( "EnergyDep" )<<"\t"
	      <<moduleHitsVector[counter][inCounter]->mother_particle().name()<<"\t"
	      <<endl;
	}
      



      //ClusteringHits(moduleHitsVector[counter], counter);
      ClusteringHits2(moduleHitsVector[counter], counter);
      //if(moduleHitsVector[counter].size() != 0){clusteringXY(moduleHitsVector[counter], counter);}
    }
}

void gdml_hit_constructor::clustering2(const std::vector<bhep::hit*>& sortedHits)
{
  /*
    Cluster the real hits (bar positions from hits) to produce hit positions.
    Also utilize the bar overlap to be able to give an even better position.
    Main jobs is done by calling clusteringXY.
  */  

  std::vector<bhep::hit*> zSortedHits = FilteringBadHits(sortedHits); 


  std::vector<bhep::hit*>::const_iterator hitIt;
  std::vector<bhep::hit*> zPlaneHits;
  std::vector<std::vector<bhep::hit*> > zPlaneHitsVector;
  //std::vector<std::vector<double> > clustered_hits;



  // Fill vectors with hits in the same module (xyxy),sorted by barPosZ.
  for (hitIt = zSortedHits.begin();hitIt != zSortedHits.end();hitIt++)
    {
      double currZ = (*hitIt)->ddata( "barPosZ" );
      double nextZ;
      rawHitsTH1F->Fill((*hitIt)->x()[2]);
      
      if(hitIt + 1 != zSortedHits.end()){ nextZ = (*(hitIt + 1))->ddata( "barPosZ" );}
      else {nextZ = currZ + 3./4. * _activeLength;}

      if(fabs(currZ-nextZ) < 3./4. * _activeLength)
	{
	  zPlaneHits.push_back((*hitIt));
	}
      else//Next is to far away
	{
	  zPlaneHits.push_back((*hitIt));

	  sort( zPlaneHits.begin(), zPlaneHits.end(), forwardSortY() );

	  zPlaneHitsVector.push_back(zPlaneHits);
	  zPlaneHits.clear();	  
	}    
    }
  /*
  for(int counter = 0; counter < zPlaneHitsVector.size(); counter++)
    {
      
      cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
	  <<"\t"<<"IsYBar"<<"\t"<<"barNum"
	  <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
	  <<"\t"<<"time"<<"\t"<<"EnergyDep"
	  <<"\t"<<"mother name"<<endl;
      for(int inCounter = 0; inCounter < zPlaneHitsVector[counter].size(); inCounter++)
	{
	  cout<<zPlaneHitsVector[counter][inCounter]->x()[0]<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->x()[1]<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->x()[2]<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->idata( "IsYBar" )<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->idata( "barNumber" )<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->ddata( "barPosZ" )<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->ddata( "barPosT" )<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->ddata( "time" )<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->ddata( "EnergyDep" )<<"\t"
	      <<zPlaneHitsVector[counter][inCounter]->mother_particle().name()<<"\t"
	      <<endl;
    }
      
    }
  */


  // Sort hits by in the same y +- 15 sorted by y.

  double tolerance = 15;

  std::vector<bhep::hit*> yPlaneHits;
  std::vector<std::vector<bhep::hit*> > yPlaneHitsVector;
  std::vector<std::vector<std::vector<bhep::hit*> > > zyPlaneHitsVector;
  std::vector<bhep::hit*>::const_iterator hitIt2;

  for(int i=0; i<zPlaneHitsVector.size();i++)
    {
      for (hitIt2 = zPlaneHitsVector[i].begin();hitIt2 != zPlaneHitsVector[i].end();hitIt2++)
	{
	  double currY = (*hitIt2)->x()[1];
	  double nextY;
	  
	  if(hitIt2 + 1 != zPlaneHitsVector[i].end()){nextY = (*(hitIt2 + 1))->x()[1];}
	  else {nextY = currY + 2*tolerance;}
	 
	  if(fabs(currY-nextY) < tolerance){yPlaneHits.push_back((*hitIt2));}
	  else//Next is to far away
	    {
	      yPlaneHits.push_back((*hitIt2));
	      sort( yPlaneHits.begin(), yPlaneHits.end(), forwardSortX() );
	      yPlaneHitsVector.push_back(yPlaneHits);
	      yPlaneHits.clear(); 
	    }    
	}
      zyPlaneHitsVector.push_back(yPlaneHitsVector);
      yPlaneHitsVector.clear();    
    }

  cout<<"after y sorting"<<endl;
  cout<<"zyPlaneHitsVector.size() "<<zyPlaneHitsVector.size()<<endl;
    
  // Sort hits by in the same y +- 15 sorted by y.

  tolerance = 15;

  std::vector<bhep::hit*> xPlaneHits;
  std::vector<std::vector<bhep::hit*> > xPlaneHitsVector;
  std::vector<std::vector<std::vector<bhep::hit*> > > yxPlaneHitsVector;
  std::vector<std::vector<std::vector<std::vector<bhep::hit*> > > > zyxPlaneHitsVector;
  //std::vector<bhep::hit*>::const_iterator hitIt3;

  for(int i=0; i<zyPlaneHitsVector.size();i++)
    {
      for(int j=0; j<zyPlaneHitsVector[i].size();j++)
	{
	  for (hitIt2 = zyPlaneHitsVector[i][j].begin();hitIt2 != zyPlaneHitsVector[i][j].end();hitIt2++)
	    {
	      double currX = (*hitIt2)->x()[0];
	      double nextX;
	      
	      if(hitIt2 + 1 != zyPlaneHitsVector[i][j].end()){nextX = (*(hitIt2 + 1))->x()[0];}
	      else {nextX = currX + 2*tolerance;}
	      
	      if(fabs(currX-nextX) < tolerance){xPlaneHits.push_back((*hitIt2));}
	      else//Next is to far away
		{
		  xPlaneHits.push_back((*hitIt2));
		  xPlaneHitsVector.push_back(xPlaneHits);
		  xPlaneHits.clear(); 
		}    
	    }
	  yxPlaneHitsVector.push_back(xPlaneHitsVector);
	  xPlaneHitsVector.clear();    
	}
      zyxPlaneHitsVector.push_back(yxPlaneHitsVector);
      yxPlaneHitsVector.clear();
    }

  for(int counter = 0; counter < zyxPlaneHitsVector.size(); counter++)
    {
      //cout<<"zyxPlaneHitsVector[counter].size() "<<zyxPlaneHitsVector[counter].size()<<endl;
      //cout<<"x"<<"\t"<<"y"<<"\t"<<"z"<<endl;
      for(int inCounter = 0; inCounter < zyxPlaneHitsVector[counter].size(); inCounter++)
	{
	  //cout<<"zyxPlaneHitsVector[counter][inCounter].size() "<<zyxPlaneHitsVector[counter][inCounter].size()<<endl;
	  for(int inInCounter = 0; inInCounter < zyxPlaneHitsVector[counter][inCounter].size(); inInCounter++)
	    {
	      
	      cout<<"zyxPlaneHitsVector[counter][inCounter][inInCounter].size() "<<zyxPlaneHitsVector[counter][inCounter][inInCounter].size()<<endl;
	      //cout<<"x"<<"\t"<<"y"<<"\t"<<"z"<<endl;
	      cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
		  <<"\t"<<"IsYBar"<<"\t"<<"barNum"
		  <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
		  <<"\t"<<"time"<<"\t"<<"EnergyDep"
		  <<"\t"<<"mother name"<<endl;
	      for(int inInInCounter = 0; inInInCounter < zyxPlaneHitsVector[counter][inCounter][inInCounter].size(); inInInCounter++)
		{
		  cout<<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->x()[0]<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->x()[1]<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->x()[2]<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->idata( "IsYBar" )<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->idata( "barNumber" )<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->ddata( "barPosZ" )<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->ddata( "barPosT" )<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->ddata( "time" )<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->ddata( "EnergyDep" )<<"\t"
		      <<((zyxPlaneHitsVector[counter])[inCounter])[inInCounter][inInInCounter]->mother_particle().name()<<"\t"
		      <<endl;

		}
	      
	      //ClusteringHits(zyxPlaneHitsVector[counter][inCounter][inInCounter], inInCounter);

      
	      //if(zyxPlaneHitsVector[counter][inCounter][inInCounter].size() != 0 &&
	      // CorrectHit(zyxPlaneHitsVector[counter][inCounter][inInCounter]))
	      if(zyxPlaneHitsVector[counter][inCounter][inInCounter].size() != 0)
		{ClusteringHits(zyxPlaneHitsVector[counter][inCounter][inInCounter], inInCounter);}
	      //{clusteringXY(zyxPlaneHitsVector[counter][inCounter][inInCounter], inInCounter);}
	    }
	  
	}
    }
  
}






std::vector<bhep::hit*> gdml_hit_constructor::FilteringBadHits(const std::vector<bhep::hit*> hits)
{

  std::vector<bhep::hit*> filteredHits;

  for(int counter = 0; counter < hits.size(); counter++)
    {
      //if(hits[counter]->ddata( "time" ) > 13.0 || hits[counter]->ddata( "EnergyDep" )< 0.1)
      if(hits[counter]->ddata( "EnergyDep" )< 0.1)
      {
        //cout<<"Removing hit from: "<<hits[inCounter]->mother_particle().name()<<endl;
        continue;
      }
      else
	{
	  filteredHits.push_back(hits[counter]);
	}

    }

  return filteredHits;
}

void gdml_hit_constructor::ClusteringHits(const std::vector<bhep::hit*> hits, int key)
{
  std::vector<bhep::hit*> filteredHits;

  //filteredHits = FilteringBadHits(hits);

  filteredHits = hits;

  //sort by time.
  sort( filteredHits.begin(), filteredHits.end(), timeSort());

  /*
  cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
      <<"\t"<<"IsYBar"<<"\t"<<"barNum"
      <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
      <<"\t"<<"time"<<"\t"<<"EnergyDep"
      <<"\t"<<"mother name"<<endl;
  for(int inCounter = 0; inCounter < filteredHits.size(); inCounter++)
    {
       cout<<filteredHits[inCounter]->x()[0]<<"\t"
	  <<filteredHits[inCounter]->x()[1]<<"\t"
	  <<filteredHits[inCounter]->x()[2]<<"\t"
	  <<filteredHits[inCounter]->idata( "IsYBar" )<<"\t"
	  <<filteredHits[inCounter]->idata( "barNumber" )<<"\t"
	  <<filteredHits[inCounter]->ddata( "barPosZ" )<<"\t"
	  <<filteredHits[inCounter]->ddata( "barPosT" )<<"\t"
	  <<filteredHits[inCounter]->ddata( "time" )<<"\t"
	  <<filteredHits[inCounter]->ddata( "EnergyDep" )<<"\t"
	  <<filteredHits[inCounter]->mother_particle().name()<<"\t"
	  <<endl;
    }
  */
    

  // Fill vectors with hits in the same time.
  double tolerance = 1; // 1 ns

  std::vector<bhep::hit*> timeHits;
  std::vector<std::vector<bhep::hit*> > timeHitsVector;

  std::vector<bhep::hit*>::const_iterator hitIt;

  for (hitIt = filteredHits.begin();hitIt != filteredHits.end();hitIt++)
    {
      double currT = (*hitIt)->ddata( "time" );
      double nextT;
      
      if(hitIt + 1 != filteredHits.end()){ nextT = (*(hitIt + 1))->ddata( "time" );}
      else {nextT = currT + 2*tolerance;}

      if(fabs(currT-nextT) < tolerance)
	{
	  timeHits.push_back((*hitIt));
	}
      else//Next is to far away
	{
	  timeHits.push_back((*hitIt));

	  if(CorrectHit(timeHits))
	    {
	      timeHitsVector.push_back(timeHits);
	    }
	  
	  timeHits.clear(); 
	}    
    }
  
  // cout<<"Hits size: "<<timeHitsVector.size()<<endl;

  cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
      <<"\t"<<"IsYBar"<<"\t"<<"barNum"
      <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
      <<"\t"<<"time"<<"\t"<<"mother name"<<endl;
  
  for(int counter = 0; counter < timeHitsVector.size(); counter++)
    {
      
      for(int inCounter = 0; inCounter < timeHitsVector[counter].size(); inCounter++)
	{
	  cout<<timeHitsVector[counter][inCounter]->x()[0]<<"\t"
	      <<timeHitsVector[counter][inCounter]->x()[1]<<"\t"
	      <<timeHitsVector[counter][inCounter]->x()[2]<<"\t"
	      <<timeHitsVector[counter][inCounter]->idata( "IsYBar" )<<"\t"
	      <<timeHitsVector[counter][inCounter]->idata( "barNumber" )<<"\t"
	      <<timeHitsVector[counter][inCounter]->ddata( "barPosZ" )<<"\t"
	      <<timeHitsVector[counter][inCounter]->ddata( "barPosT" )<<"\t"
	      <<timeHitsVector[counter][inCounter]->ddata( "time" )<<"\t"
	      <<timeHitsVector[counter][inCounter]->mother_particle().name()<<"\t"
	      <<endl;
	}
      
      
      if(timeHitsVector[counter].size() != 0){clusteringXY(timeHitsVector[counter], counter);}
    }
  
}

void gdml_hit_constructor::ClusteringHits2(const std::vector<bhep::hit*> hits, int key)
{
  std::vector<bhep::hit*> filteredHits;

  filteredHits = FilteringBadHits(hits);
  //filteredHits = hits;

  //sort by time.
  sort( filteredHits.begin(), filteredHits.end(), timeSort());

  // Fill vectors with hits in the same time.
  double tolerance = 1; // 1 ns

  std::vector<bhep::hit*> timeHits;
  std::vector<std::vector<bhep::hit*> > timeHitsVector;

  std::vector<bhep::hit*>::const_iterator hitIt;

  for (hitIt = filteredHits.begin();hitIt != filteredHits.end();hitIt++)
    {
      double currT = (*hitIt)->ddata( "time" );
      double nextT;
      
      if(hitIt + 1 != filteredHits.end()){ nextT = (*(hitIt + 1))->ddata( "time" );}
      else {nextT = currT + 2*tolerance;}

      if(fabs(currT-nextT) < tolerance)
	{
	  timeHits.push_back((*hitIt));
	}
      else//Next is to far away
	{
	  timeHits.push_back((*hitIt));

	  // if(CorrectHit(timeHits))
	  //{
	  timeHitsVector.push_back(timeHits);
	  //}
	  
	  timeHits.clear(); 
	}    
    }

  // Now we have a vector with hits in the same z-plane, close in time. 
  // From this, create hits and "ghosthits" if more than 1 x or y hit exists.

  std::vector<bhep::hit*> X, Y;
  std::vector<std::vector<bhep::hit*> > xHitsVector, yHitsVector;
  
  // Sort out the x and y bar hits.
  for(int counter = 0; counter < timeHitsVector.size(); counter++)
    {
      for(int inCounter = 0; inCounter < timeHitsVector[counter].size(); inCounter++)
	{
	  if( timeHitsVector[counter][inCounter]->idata( "IsYBar" ) == 0){X.push_back(timeHitsVector[counter][inCounter]);}
	  else {Y.push_back(timeHitsVector[counter][inCounter]);}
	}
      xHitsVector = NextCluster(X, 2, "barNumber");//, xHitsVector);
      yHitsVector = NextCluster(Y, 2, "barNumber");//, yHitsVector);

      std::vector<std::vector<bhep::hit*> > totHitsVector;
      std::vector<bhep::hit*> totHits;
      
      //combine all x with all y.
      
      for(int i=0;i<xHitsVector.size();i++)
	{
	  for(int j=0;j<yHitsVector.size();j++)
	    {
	      totHits = xHitsVector[i];
	      totHits.insert(totHits.end(),yHitsVector[j].begin(),yHitsVector[j].end());
	      totHitsVector.push_back(totHits);
	      totHits.clear();
	    }
	}

      for(int counter=0; counter<totHitsVector.size();counter++)
	{
	  cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
	      <<"\t"<<"IsYBar"<<"\t"<<"barNum"
	      <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
	      <<"\t"<<"time"<<"\t"<<"mother name"<<endl;
	  cout<<"totHitsVector[counter].size() "<<totHitsVector[counter].size()<<endl;
	  for(int inCounter=0;inCounter<totHitsVector[counter].size();inCounter++)
	    {
	  cout<<totHitsVector[counter][inCounter]->x()[0]<<"\t"
	      <<totHitsVector[counter][inCounter]->x()[1]<<"\t"
	      <<totHitsVector[counter][inCounter]->x()[2]<<"\t"
	      <<totHitsVector[counter][inCounter]->idata( "IsYBar" )<<"\t"
	      <<totHitsVector[counter][inCounter]->idata( "barNumber" )<<"\t"
	      <<totHitsVector[counter][inCounter]->ddata( "barPosZ" )<<"\t"
	      <<totHitsVector[counter][inCounter]->ddata( "barPosT" )<<"\t"
	      <<totHitsVector[counter][inCounter]->ddata( "time" )<<"\t"
	      <<totHitsVector[counter][inCounter]->mother_particle().name()<<"\t"
	      <<endl;
	    }


	  if(totHitsVector[counter].size() != 0 && CorrectHit(totHitsVector[counter]))
	    {clusteringXY(totHitsVector[counter], counter);}
	}
      
    }

}
//void gdml_hit_constructor::NextCluster(std::vector<bhep::hit*> inHits,
std::vector<std::vector<bhep::hit*> > gdml_hit_constructor::NextCluster(std::vector<bhep::hit*> inHits,
									double tolerance, 
									string data)
//std::vector<std::vector<bhep::hit*> > hitsVector)
{
  std::vector<bhep::hit*>::const_iterator hitIt;
  std::vector<bhep::hit*> hits;
  std::vector<std::vector<bhep::hit*> > hitsVector;

  for (hitIt = inHits.begin();hitIt != inHits.end();hitIt++)
    {
      //cout<<"In for"<<endl;
      double curr = (*hitIt)->idata( "barNumber" );
      double next;
      
      if(hitIt + 1 != inHits.end()){ next = (*(hitIt + 1))->idata( "barNumber" );}
      else {next = curr + 2*tolerance;}
      
      if(fabs(curr-next) < tolerance){hits.push_back((*hitIt));}
      else//Next is to far away
	{
	  hits.push_back((*hitIt));
	  hitsVector.push_back(hits);
	  //cout<<"Filled"<<endl;
	  hits.clear(); 
	}    
    } 

  return hitsVector;
} 



bool  gdml_hit_constructor::CorrectHit(const std::vector<bhep::hit*> hits)
{
  /*
    Ensure that the vector containes atleast one x and one y bar hit.
  */

  int x = 0;
  int y = 0;
  bool correct = false;

  for(int counter = 0; counter < hits.size(); counter++)
    {
      if(hits[counter]->idata( "IsYBar" ) ==0)
	{
	  x++;
	}
      else
	{
	  y++;
	}
    }

  if(x > 0 && y >0)
    {
      correct = true;
    }
  else
    {
      cout<<"Incorrect hit"<<endl;
      cout<<x<<endl;
      cout<<y<<endl;
      cout<<hits.size()<<endl;
    }
  
  return correct;
}


void gdml_hit_constructor::clusteringXY(const std::vector<bhep::hit*> hits, int key)
{
  /*
    Cluster the real hits (bar positions from hits) to produce hit positions.
    Also utilize the bar overlap to be able to give an even better position.
    Results in _voxel being filled.
  */  
  std::vector<bhep::hit*> X, Y;

  std::vector<bhep::hit*> filteredHits;

  double z = 0;
  /*
  cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
      <<"\t"<<"IsYBar"<<"\t"<<"barNum"
      <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
      <<"\t"<<"time"<<"\t"<<"EnergyDep"
      <<"\t"<<"mother name"<<endl;
  */
  
  for(int inCounter = 0; inCounter < hits.size(); inCounter++)
    {
      /*
      cout<<hits[inCounter]->x()[0]<<"\t"
	  <<hits[inCounter]->x()[1]<<"\t"
	  <<hits[inCounter]->x()[2]<<"\t"
	  <<hits[inCounter]->idata( "IsYBar" )<<"\t"
	  <<hits[inCounter]->idata( "barNumber" )<<"\t"
	  <<hits[inCounter]->ddata( "barPosZ" )<<"\t"
	  <<hits[inCounter]->ddata( "barPosT" )<<"\t"
	  <<hits[inCounter]->ddata( "time" )<<"\t"
	  <<hits[inCounter]->ddata( "EnergyDep" )<<"\t"
	  <<hits[inCounter]->mother_particle().name()<<"\t"
	  <<endl;
      */
      
      /*
      if(hits[inCounter]->ddata( "time" ) > 0.1 || hits[inCounter]->ddata( "EnergyDep" )< 0.1)
      {
        cout<<"Removing hit from: "<<hits[inCounter]->mother_particle().name()<<endl;
        continue;
      }
      */
      filteredHits.push_back(hits[inCounter]);

      z+=hits[inCounter]->ddata( "barPosZ" );
      
      if( hits[inCounter]->idata( "IsYBar" ) == 0){X.push_back(hits[inCounter]);}
      else {Y.push_back(hits[inCounter]);}
    }
  
  z= z/filteredHits.size();
  int vox_x = -1;
  int vox_y = -1;

  //cout<<"More than 4? "<<X.size()<<" "<<Y.size()<<endl;
  
  if(X.size() != 0)
    {
      vox_x = calculate_vox_no(X);
 
    }
  if(Y.size() != 0)
    {
      vox_y = calculate_vox_no(Y);
 
    }

  cout<<"vox_x "<<vox_x<<endl;
  cout<<"vox_y "<<vox_y<<endl;

  int vox_num = vox_x + vox_y*_nVoxX;
  cout<<"vox_num "<<vox_num<<endl;
  cout<<"z "<<z<<endl;


  //for the whole vector.

  for(int cnt = 0; cnt<filteredHits.size();cnt++)
    { 
      if ( vox_num >= 0){
	_voxels[z].insert( pair<int,bhep::hit*>(vox_num,filteredHits[cnt]) );
	clusteredHitsTH1F->Fill(filteredHits[cnt]->x()[2]);
      }
    }
    
}

int gdml_hit_constructor::calculate_vox_no(std::vector<bhep::hit*> hits)
{
  /*
    Take all the hits in a module (4z planes) and calculate the voxel number from it.
    Takes in only X or only Y plane hits. Does not yet handle multiple hits per plane.

  */

  int vox_num = -1;
  int currBarNum = hits[0]->idata( "barNumber" );

  if(hits.size() == 1)
    {
      //int currBarNum = hits[0]->idata( "barNumber" );

      if(currBarNum % 2 == 0) // If barNum even then back bar
	{
	  vox_num = 2*currBarNum;
	}
      else{vox_num = 2*(currBarNum -1)+2 ;}
    }
  else
    {
      int frontBarNum = -1;
      int backBarNum = -1;

      for(int cnt = 0; cnt<hits.size(); cnt++)
	{
	  currBarNum = hits[cnt]->idata( "barNumber" );
	  if(currBarNum % 2 == 0){backBarNum = currBarNum/2;}
	  else{frontBarNum = (int) (currBarNum -1)/2;}
	}
      vox_num = 2*backBarNum + 2*frontBarNum + 1;
    }

  // cout<<"vox_num: "<<vox_num<<endl;

 return vox_num;
}



void gdml_hit_constructor::construct_hits(std::vector<bhep::hit*>& rec_hit)
{
  //takes the voxels which have been filled and make
  //rec_hit objects out of them.

  std::map<double,std::multimap<int,bhep::hit*> >::iterator vIt;
  std::multimap<int,bhep::hit*>::iterator vIt2;

  for (vIt = _voxels.begin();vIt != _voxels.end();vIt++){

    while ( vIt->second.size() != 0 ){
      //set second iterator to first filled voxel in layer.
      vIt2 = vIt->second.begin();

      bhep::hit* vhit = get_vhit( vIt2->first, vIt->first, vIt->second );
      
      if ( vhit != NULL ){
	rec_hit.push_back( vhit );
	
	//std::cout<<"In construct_hits in hit_constructor "<<"x = "<<vhit->x()[0]<<",y = "<<vhit->x()[1]<<", z = "<<vhit->x()[2]<<std::endl;
      }
      vIt->second.erase( vIt2->first );
    }
  }
  
}

bhep::hit* gdml_hit_constructor::get_vhit(int vox, double z,
				     const std::multimap<int,bhep::hit*>& map1)
{
  //Makes a rec_hit from the voxel position and adds the relevant points.
  bhep::hit* returnPointer;  

  double totEng = 0., muProp = 0.; //these will be done in rec_hit eventually.
  vdouble X, Y, Z, E, T; //annoying but again all in rec_hit class.
  double proptime=9999999.9, vlight = 299792458. / 1.6;
  double meanvoxtime;

  double smearingFactor = 1;
  double dt, dtx, dty;
  double xedge = _detectorX/2.;
  double yedge = _detectorY/2.;
  double xE1, xE2, yE1, yE2;
  double Xphot, Yphot;

  vector<double> barPosX;
  vector<double> barPosY;
  double sumBarPosX = 0;
  double sumBarPosY = 0;
  double barX;
  double barY;

  vstring mother_particle;

  double propSpeed = 170000000.0;
  double timeSmearingFactor = 0.01;


  bhep::hit* vhit = new bhep::hit( "tracking" );

  vhit->add_property( "voxel", vox );


  
  // Pushback all the relevant informatiomn to the voxels. 
  // Also calculate the average position of the hits in each voxel.
  std::multimap<int,bhep::hit*>::const_iterator hIt;


  /*
  cout<<"x"<<"\t"<<"y"<<"\t"<<"z"
      <<"\t"<<"IsYBar"<<"\t"<<"barNum"
      <<"\t"<<"barPosZ"<<"\t"<<"barPosT"
      <<"\t"<<"time"<<"\t"<<"EnergyDep"
      <<"\t"<<"mother name"<<endl;
  */
  for (hIt = map1.equal_range(vox).first;hIt != map1.equal_range(vox).second;hIt++)
    {  
      /*
      cout<<(*hIt).second->x()[0]<<"\t"
	  <<(*hIt).second->x()[1]<<"\t"
	  <<(*hIt).second->x()[2]<<"\t"
	  <<(*hIt).second->idata( "IsYBar" )<<"\t"
	  <<(*hIt).second->idata( "barNumber" )<<"\t"
	  <<(*hIt).second->ddata( "barPosZ" )<<"\t"
	  <<(*hIt).second->ddata( "barPosT" )<<"\t"
	  <<(*hIt).second->ddata( "time" )<<"\t"
	  <<(*hIt).second->ddata( "EnergyDep" )<<"\t"
	  <<(*hIt).second->mother_particle().name()<<"\t"
	  <<endl;
      */

      //get propagation time.
      double xt1 = (xedge+(*hIt).second->x()[0])/propSpeed;
      double xt2 = (xedge-(*hIt).second->x()[0])/propSpeed;

      double yt1 = (yedge+(*hIt).second->x()[1])/propSpeed;
      double yt2 = (yedge-(*hIt).second->x()[1])/propSpeed;

      // smear time
      xt1 = xt1 + _ranGen.Gaus(0,timeSmearingFactor * xt1);
      xt2 = xt2 + _ranGen.Gaus(0,timeSmearingFactor * xt2);

      yt1 = yt1 + _ranGen.Gaus(0,timeSmearingFactor * yt1);
      yt2 = yt2 + _ranGen.Gaus(0,timeSmearingFactor * yt2);

      // Convert back to position.

      double x = (xt1-xt2)/2.0 * propSpeed;
      double y = (yt1-yt2)/2.0 * propSpeed;

      if((*hIt).second->idata( "IsYBar" ) == 0)
	{
	  barPosX.push_back((*hIt).second->ddata( "barPosT" ));
	  sumBarPosX+=(*hIt).second->ddata( "barPosT" );
	  /*
	  cout<<"x-bar x="<<(*hIt).second->ddata( "barPosT" )<<endl;  
	  cout<<"x-bar y="<<y<<endl; 
 	  cout<<"x-bar real y="<<(*hIt).second->x()[1]<<endl;
	  */
	  barPosY.push_back(y);
	  sumBarPosY+=y;
	}
      else
	{
	  barPosY.push_back((*hIt).second->ddata( "barPosT" ));
	  sumBarPosY+=(*hIt).second->ddata( "barPosT" );
	  /*
	  cout<<"y-bar y="<<(*hIt).second->ddata( "barPosT" )<<endl;  
	  cout<<"y-bar x="<<x<<endl; 
	  cout<<"y-bar real x="<<(*hIt).second->x()[0]<<endl;
	  */
	  barPosX.push_back(x);
	  sumBarPosX+=x;
	}
      mother_particle.push_back((*hIt).second->mother_particle().name());
      X.push_back( (*hIt).second->x()[0] );
      Y.push_back( (*hIt).second->x()[1] );
      Z.push_back( (*hIt).second->x()[2] );
      digitizedHitsTH1F->Fill((*hIt).second->x()[2]);
      T.push_back( (*hIt).second->ddata( "time" ) );
      E.push_back( (*hIt).second->ddata( "EnergyDep" ) );
      totEng += (*hIt).second->ddata( "EnergyDep" );
      proptime = (*hIt).second->ddata( "time" )  < proptime ?  
	(*hIt).second->ddata( "time" )  : proptime;
      if ( (*hIt).second->mother_particle().name() == "mu+" ||
	   (*hIt).second->mother_particle().name() == "mu-" ){
	if ( (*hIt).second->mother_particle().fetch_sproperty("CreatorProcess")=="none" )
	  muProp++;
      } else if ( (*hIt).second->mother_particle().name() == "lepton_shower" )
	muProp += 0.5;
      
      //std::cout<<(*hIt).second->x()[0]<<"\t"<<(*hIt).second->x()[1]<<"\t"<<(*hIt).second->x()[2]<<"\t"
      //	       <<(*hIt).second->ddata( "EnergyDep" )<<std::endl;
    }
  barX=sumBarPosX/barPosX.size();
  barY=sumBarPosY/barPosY.size();

  Point3D hitPos(barX , barY, z );
  vhit->set_point( hitPos );

  //Do attenuations.

  //Assume equal amounts of energy from both views and
  // equal energy flow in both directions along strip.
  //cout<<"totEng: "<<totEng<<endl;
  //xE1 = xE2 = yE1 = yE2 = totEng/2;

  xE1 = xE2 = yE1 = yE2 = totEng/4;

  xeTH1F->Fill(xE1+xE2);
  yeTH1F->Fill(yE1+yE2);

  //cout<<"attLength: "<<_attLength<<endl;
  /*
  xE1 = xE1 * exp(-(xedge + fabs(barX))/_attLength);
  xE2 = xE2 * exp(-(3*xedge-fabs(barX))/_attLength);
  yE1 = yE1 * exp(-(yedge + fabs(barY))/_attLength);
  yE2 = yE2 * exp(-(3*yedge-fabs(barY))/_attLength);
  */
  xE1 = xE1 * exp(-(xedge - fabs(barX))/_attLength);
  xE2 = xE2 * exp(-(xedge + fabs(barX))/_attLength);
  yE1 = yE1 * exp(-(yedge - fabs(barY))/_attLength);
  yE2 = yE2 * exp(-(yedge + fabs(barY))/_attLength);


  xeAttTH1F->Fill(xE1+xE2);
  yeAttTH1F->Fill(yE1+yE2);
  /*
  dtx = (xedge + fabs(barX)) < (3*xedge - fabs(barX)) ?
    (xedge + fabs(barX))/vlight : (3*xedge - fabs(barX))/vlight;
  dty = (yedge + fabs(barY)) < (3*yedge - fabs(barY)) ?
    (yedge + fabs(barY))/vlight : (3*yedge - fabs(barY))/vlight;
  */
  dtx = (xedge - fabs(barX)) < (xedge + fabs(barX)) ?
    (xedge - fabs(barX))/vlight : (xedge + fabs(barX))/vlight;
  dty = (yedge - fabs(barY)) < (yedge + fabs(barY)) ?
    (yedge - fabs(barY))/vlight : (yedge + fabs(barY))/vlight;


  dt = dtx < dty ? dtx : dty;

  //double xE = xE1 + _ranGen.Gaus( 0, smearingFactor * xE1 )
  //  + xE2 + _ranGen.Gaus( 0, smearingFactor * xE2 );
  double xE = xE1 + _ranGen.PoissonD(smearingFactor * xE1)
    + xE2 + _ranGen.PoissonD(smearingFactor * xE2);
  //double yE = yE1 + _ranGen.Gaus( 0, smearingFactor * yE1 )
  //  + yE2 + _ranGen.Gaus( 0, smearingFactor * yE2 );
  double yE = yE1 + _ranGen.PoissonD(smearingFactor * yE1)
    + yE2 + _ranGen.PoissonD(smearingFactor * yE2);

  xeSmearTH1F->Fill(xE);
  yeSmearTH1F->Fill(yE);

  //if ( fabs(z) > (_detectorLength + _vertexDetdepth)/2. && 
  //     xE < _minEng && yE < _minEng )
  if ( fabs(z) > (_detectorLength)/2. && xE < _minEng && yE < _minEng )
    {
      delete vhit;
      returnPointer = NULL;
    }
  else
    {
      totEng = xE + yE;
      proptime += dt;
      //cout<<"whit properties added"<<endl;
      vhit->add_property( "mother_particle", mother_particle );
      
      vhit->add_property( "TotalEng", totEng );
      vhit->add_property( "XEng", xE );
      vhit->add_property( "YEng", yE );
      vhit->add_property( "MuonProportion", muProp );
      vhit->add_property( "NoPoints", (int)X.size() );
      vhit->add_property( "Xpoint", X );
      vhit->add_property( "Ypoint", Y );
      vhit->add_property( "Zpoint", Z );
      vhit->add_property( "Epoint", E );
      vhit->add_property( "HitTime", proptime);
      
      returnPointer = vhit;
    }
  
  return returnPointer;
}
				     

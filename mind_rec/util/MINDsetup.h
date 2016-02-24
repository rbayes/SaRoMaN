/* -*- mode: c++ -*- */
#ifndef _mind_setup___
#define _mind_setup___

#include <recpack/RecPackManager.h>
#include <mind/SetupSk.h>
#include <mind/MINDfieldMapReader.h>
#include <mind/DeDxMap.h>

//#include <string.h>

//USING THE GDML PARSED FILE
#include <map>
#include <iostream>
#include <fstream>

using namespace Recpack;

class MINDsetup: public NSetupSk{

public:
    
  MINDsetup();
  
  virtual ~MINDsetup();
    
  Setup& setup();
   
  // void init(bhep::gstore store,bhep::sstore,
  // 	    bhep::prlevel level=bhep::NORMAL);
  void init(bhep::gstore store, bhep::prlevel level=bhep::NORMAL);
  
  //info to build virtual planes
  
  double getPlaneX(){return MIND_x;};
  double getPlaneY(){return MIND_y;};
  double getPlaneZ(){return MIND_z;};
  EVector getXaxis(){return xaxis;};
  EVector getYaxis(){return yaxis;};
  EVector getZaxis(){return zaxis;};
  string getMeasType(){return meastype;};
  int getMeasDim(){return meas_dim;} 
  EVector getResolution(){return resolution;} 
  EMatrix getCov(){return cov;};
  double get_Fe_prop(){return _wFe;}
  double& getDeDx(){return de_dx;}
  void setDeDx(double d){de_dx = d;}
  //EVector getBField(EVector pos){return BFieldMap.vector(pos);}
  EVector getBField(EVector pos);

  EVector getRawBField(EVector pos);

  double getPieceWidth() {return _pieceWidth;}

  std::map<dict::Key,vector<double> > getModuleDataMap() {return _moduleDataMap;}

  // vol_name, module position z, module size z, wFe, MagFieldScale
  std::map<dict::Key,vector<double> > _moduleDataMap;

  Volume* mother;
  Volume* det;
  vector<Volume*> detVector;
  Volume* vdet;
  Volume* Fe_slab;

  // z max and min in the detector volume.
  double getZMax() {return _detector_z_max;}
  double getZMin() {return _detector_z_min;}
  double getBScaleAvr() {return _detector_Bscale_avr;}

  double _detector_z_max;
  double _detector_z_min;
  double _detector_Bscale_avr;

 
protected:
    
  void readParam();
  void createGeom();
  void add_slab(int plane, const dict::Key det_vol);
  // void add_slab(int plane, Volume& det_vol);
  void setResolution();
  void addProperties();
  

protected:
  
  // parameter store

  bhep::gstore _pstore;
  
  // detector axis

  EVector xaxis,yaxis,zaxis;
    
  // -------------------------------------------------------------//
  //                       |  DIMENSIONS |                        //
  // -------------------------------------------------------------//
    
  //--------------------------- VOLUMES -------------------------//
     
  int    OctGeom;
  double MOTHER_x, MIND_x, VERT_x;
  double MOTHER_y, MIND_y, VERT_y;
  double MOTHER_z, MIND_z, VERT_z;
  double MOTHER_earh, EAR_height;
  double MOTHER_earw, EAR_width;
  double IRON_z, SCINT_z, AIR_z, Al_z;
  double rel_denAS, rel_denSI;//AIR/Scint, Scint/Fe.
  int nScint;
  int _npieces;
  double _pieceWidth;

  // -------------------------------------------------------------//
  //                         |  PHYSICS |                         //
  // -------------------------------------------------------------//

  //------------------------ MAGNETIC FIELD ----------------------//
    
  double B_int;
  EVector BField;
  //double _fieldScale;
  //MINDfieldMapReader BFieldMap;
  vector<MINDfieldMapReader*> BFieldMapVec;
  MINDfieldMapReader _generalBFieldMap;

  //-------------------------------------------------------------//
  
  //------------------- PROPERTIES OF MATERIALS -----------------//
    
  double X0Fe, X0Sc, X0AIR;//members for if/when geom more strict.
  double de_dx;
  double _wFe;

  vector<double> X0EffVec;

  //for de/dx map
  double de_dx_scint;
  double de_dx_fe;
  //double de_dx_min;
  vector<DeDxMap*> de_dx_map_vec;
  //DeDxMap* _de_dx_map;
  //DeDxMap* _de_dx_map_scint;
  std::string Bmap;
  
  EVector _zaxis;
  
  //-------------------------------------------------------------//
  

  //------------------------- MEASUREMENTS ----------------------//
     
  int meas_dim;
  string meastype;
  EVector resolution;
  EMatrix cov;
  double resx,resy,resz;
  double StepSize;



  //------------------------- GDML PARSED FILE ----------------------//

  string _gdml_parsed_path;

  std::map<string,std::vector<double> > _gdml_solid_map;

  std::map<string,std::vector<double> > _gdml_pos_map;

};



#endif 

#include <mind/MINDsetup.h>
#include <recpack/string_tools.h>
#include <recpack/dictionary.h>
#include <string>
#include <algorithm>
#include <mind/MINDplate.h>
#include <mind/EMINDplate.h>
#include <mind/DeDxMap.h>

#include <bhep/system_of_units.h>



//*************************************************************
MINDsetup::MINDsetup() {
//*************************************************************
  
}


//*************************************************************
MINDsetup::~MINDsetup() {
//*************************************************************

}


//*************************************************************
void MINDsetup::init(bhep::gstore pstore, bhep::prlevel level) {
//*************************************************************
    
  _msetup=bhep::messenger(level);
  
  _msetup.message("++MINDsetup Messenger generated++",bhep::VERBOSE);
  
  _pstore=pstore;
  
  readParam();

  //--------------- generate recpack setup -----------//
    
  _gsetup.set_name("main");
  
  // create volumes and surfaces
  
  createGeom();
  
  // define detector resolutions

  setResolution();
  
  // add properties to volumes and surfaces
  
  addProperties();

  // std::cout << _gsetup.volume("Detector").parameter("BField") << std::endl;
  // dict::mixdictionary smell = _gsetup.volume_properties("IRON_plane0");
  
  // const double XX = _gsetup.volume_properties("IRON_plane0").retrieve(thing_name);
  // std::cout << _gsetup.volume_properties("IRON_plane0").retrieve(thing_name) << std::endl;
  
  //     std::cout << _gsetup << std::endl;
  
  _msetup.message("++ Setup has been generated !! ++",bhep::VERBOSE);
  
  //_msetup.message("MIND Setup:", _gsetup,bhep::VERBOSE);
}


//*************************************************************
Setup& MINDsetup::setup() {
//*************************************************************
     
  return _gsetup;

}

//*************************************************************
void MINDsetup::createGeom(){
//*************************************************************    
    
  _msetup.message("+++ CreatGeom function +++",bhep::VERBOSE);
  
  //----- axes for definition of volumes and surfaces ----//
  
  xaxis=EVector(3,0); xaxis[0] = 1.; 
  yaxis=EVector(3,0); yaxis[1] = 1.; 
  zaxis=EVector(3,0); zaxis[2] = 1.;
  
  //----------- Create a mandatory mother volume ---------//
  
  EVector pos(3,0); pos[2]= VERT_z/2.;
  EVector vpos(3,0); vpos[2] = -MIND_z/2.;
  /*
    Volume* mother = new Box(pos,xaxis,yaxis,
    MOTHER_x/2,MOTHER_y/2,MOTHER_z/2);
  */
  Volume* mother;
  if(OctGeom==1)
    mother = new MINDplate(pos,xaxis,yaxis,
			   MOTHER_x/2,MOTHER_y/2,MOTHER_z/2,
			   MOTHER_earw,MOTHER_earh);
  else if(OctGeom==2)
    mother = new EMINDplate(pos,xaxis,yaxis,
			    MOTHER_x/2,MOTHER_y/2,MOTHER_z/2,
			    MOTHER_earw,MOTHER_earh);

  else if(OctGeom==3)
    mother = new Box(pos, zaxis, xaxis, MOTHER_z/2., MOTHER_x/2, MOTHER_y/2);
 
  else
    mother = new Tube(pos, zaxis,MOTHER_z/2, MOTHER_x/2); 

  _msetup.message("Mother volume generated",bhep::VERBOSE);

  // add mother volume
  
  _gsetup.set_mother(mother);
   
  _msetup.message("Mother added to setup",bhep::VERBOSE);

  // Create detector volume

  const dict::Key vol_name = "Detector";
  const dict::Key vert_name = "VertDetector";
  // Volume* det = new Box(pos,xaxis,yaxis,MIND_x/2,MIND_y/2,MIND_z/2);

  Volume* det;
  Volume* vdet;
  if(OctGeom==1)
    det = new MINDplate(pos,xaxis,yaxis,MIND_x/2,MIND_y/2,MIND_z/2,
				EAR_width, EAR_height);
  else if(OctGeom==2)
    det = new EMINDplate(pos,xaxis,yaxis,
			    MIND_x/2,MIND_y/2,MIND_z/2,
			    EAR_width,EAR_height);
  else if(OctGeom==3)
    det = new Box(pos, zaxis, xaxis, MOTHER_z/2., MOTHER_x/2, MOTHER_y/2);
  else{
    det = new Tube(pos,zaxis,MIND_z/2,MIND_x/2);
    if(VERT_z > 0)
      vdet = new Box(vpos, zaxis, xaxis, VERT_z/2., VERT_x/2., VERT_y/2.);
  }
  _msetup.message("MIND volume generated",bhep::VERBOSE);
 
  // add volume
  
  _gsetup.add_volume("mother",vol_name,det);
  if(VERT_z > 0)
    _gsetup.add_volume("mother",vert_name,vdet);
  // _gsetup.set_volume_property(vol_name,"X0",X0AIR);
  //Introduce IRON scintillator sandwiches.
  // int nplanes = (int)( MIND_z / (IRON_z + nScint * SCINT_z) );
  
  // for (int iplane = 0;iplane < _npieces;iplane++) {
  
  //     add_slab(iplane, vol_name);
  
  //   }
  
  
  
}

//*************************************************************
void MINDsetup::add_slab(int plane, const dict::Key det_vol){
  //*************************************************************
  
  // EVector plane_pos(3,0);
  //Names for particular sections
  
  //Define positions
  // double slab_width = IRON_z + nScint * SCINT_z;
  double mind_front = -MIND_z/2;
  
  //IRON
  EVector fe_pos(3,0);
  const dict::Key iron_name = "IRON_plane"+bhep::to_string(plane);

  fe_pos[2] = mind_front + plane*_pieceWidth + IRON_z/2;
  
  // Volume* Fe_slab = new Box(fe_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,IRON_z/2);
  
  Volume* Fe_slab = new MINDplate(fe_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,
				  IRON_z/2, EAR_width, EAR_height);
  
  _gsetup.add_volume(det_vol,iron_name,Fe_slab);
  
  _gsetup.set_volume_property(iron_name,"X0",X0Fe);
  
  //SCINT
  EVector scint_pos(3,0);
  Volume *Sc_slab[nScint];
  for (int iscint = 0;iscint < nScint;iscint++){
    
    const dict::Key scint_name =
      "SCINT_plane"+bhep::to_string(plane)+"_"+bhep::to_string(iscint);

    scint_pos[2] = mind_front + plane*_pieceWidth + IRON_z
      + (iscint+1)*AIR_z + SCINT_z/2 + Al_z * iscint;
    
    // Sc_slab[iscint] = new Box(scint_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,SCINT_z/2);
    
    Sc_slab[iscint] = new MINDplate(scint_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,
				  SCINT_z/2, EAR_width, EAR_height);

    _gsetup.add_volume(det_vol,scint_name,Sc_slab[iscint]);

    _gsetup.set_volume_property(scint_name,"X0",X0Sc);

  }

  //Scintillator.
  // plane_pos[2] = mind_front + _pieceWidth * plane + IRON_z + SCINT_z/2;
//   const dict::Key scint_name = "SCINT_plane"+bhep::to_string(plane_pos[2]);

//   Volume* Sc_slab = new Box(plane_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,SCINT_z/2);

//   _gsetup.add_volume(det_vol,scint_name,Sc_slab);

//   _gsetup.set_volume_property(scint_name,"X0",X0Sc);
  
//   //IRON
//   plane_pos[2] = mind_front + slab_width * plane + IRON_z/2;

//   Volume* Fe_slab = new Box(plane_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,IRON_z/2);

//   _gsetup.add_volume(det_vol,iron_name,Fe_slab);

//   _gsetup.set_volume_property(iron_name,"X0",X0Fe);
 
}

//*************************************************************
void MINDsetup::setResolution(){
//*************************************************************    
        
    /*
      
    */
  
  _msetup.message("+++ setResolution function +++",bhep::VERBOSE);

  resolution = EVector(meas_dim,0);
  resolution[0] = resx*mm;  
  resolution[1] = resy*mm;  
  
    
  // cov of measurements

  cov = EMatrix(meas_dim,meas_dim,0);
  for (size_t i = 0; i < (size_t)meas_dim; i++) 
    cov[i][i] = resolution[i]*resolution[i];
    
}



//*************************************************************
void MINDsetup::addProperties(){
//*************************************************************    
  
  _msetup.message("+++ addProperties function +++",bhep::VERBOSE);

  //-------------------- magnetic field ------------------//
    
  BField = EVector(3,0);
  BField[1] = B_int;

  _zaxis = EVector(3,0);
  _zaxis[2]=1;
  
  // _gsetup.set_volume_property("mother","BField",BField);
  const dict::Key vol_name = "Detector";
  const dict::Key vert_name = "VertDetector";
  //  _msetup.message("+++B Field added to MOTHER:",BField,bhep::VERBOSE);
  // step = 1*cm;
  if(OctGeom==3)
    _gsetup.set_volume_property_to_sons("mother","BField",BField);
  else { 
    _gsetup.set_volume_property_to_sons(vol_name,RP::BFieldMap,BFieldMap);
    if(VERT_z > 0.0)
      _gsetup.set_volume_property_to_sons(vert_name,"BField",BField);
  }
  // _gsetup.set_volume_property_to_sons("mother","de_dx",de_dx);

  //Instead of fixed de_dx, the energy deposition ditribution map 
   
  _de_dx_map = new DeDxMap(de_dx_min*MeV/mm);
  _de_dx_map_scint = new DeDxMap(de_dx_scint*MeV/mm);
  _gsetup.set_volume_property_to_sons(vol_name,RP::de_dx_map,*_de_dx_map);
  if(VERT_z > 0.0)
    _gsetup.set_volume_property_to_sons(vert_name,RP::de_dx_map,*_de_dx_map_scint);
  _gsetup.set_volume_property_to_sons("mother",RP::SurfNormal,_zaxis);
  // _gsetup.set_volume_property_to_sons("mother",RP::StepSize,step);
  // _gsetup.set_volume_property_to_sons(vol_name,"BField",BField);
//   _gsetup.set_volume_property_to_sons(vol_name,"de_dx",de_dx);
  // _gsetup.set_volume_property_to_sons(vol_name,RP::SurfNormal,_zaxis);

  // const dict::Key vol_name = "Detector";
  // _gsetup.set_volume_property("mother","X0",X0AIR);
  // _gsetup.set_volume_property(vol_name,"X0",X0AIR);
    
  // const dict::Key vol_name = "Detector";
  // _gsetup.set_volume_property(vol_name,"BField",BField);
  // _gsetup.set_volume_property(vol_name,"BFieldMap",BFieldMap);
  _msetup.message("+++B Field added to MIND:","BFieldMap",bhep::VERBOSE);
  
  if(StepSize){
    _gsetup.set_volume_property(vol_name,"StepSize",StepSize);
    if(VERT_z > 0.0)
      _gsetup.set_volume_property(vert_name,"StepSize",StepSize);
  }
  _gsetup.set_volume_property(vol_name,"X0",X0Eff);
  if(VERT_z > 0.0)
    _gsetup.set_volume_property(vert_name,"X0",X0Sc);
//   _msetup.message("+++X0 added to MIND:",X0,bhep::VERBOSE);

//   // _gsetup.set_volume_property(vol_name,"de_dx",de_dx);
//   _msetup.message("+++de/dx added to MIND:",de_dx,bhep::VERBOSE);
  
 
}





void MINDsetup::readParam(){

    
    // -------------------------------------------------------------//
    //                       |  DIMENSIONS |                        //
    // -------------------------------------------------------------//
    
    bhep::prlevel c = bhep::VERBOSE;
    
    OctGeom = _pstore.find_istore("IsOctagonal") ? 
      _pstore.fetch_istore("IsOctagonal"): 1;

    MIND_x = _pstore.fetch_dstore("MIND_x") * m; 
      
    _msetup.message("MIND height:",MIND_x/cm,"cm",c);
    
    MIND_y =  _pstore.fetch_dstore("MIND_y") * m;
      
    _msetup.message("MIND width:",MIND_y/cm,"cm",c);
    
    MIND_z =  _pstore.fetch_dstore("MIND_z") * m;
      
    _msetup.message("MIND length:",MIND_z/cm,"cm",c);
    
    EAR_height = _pstore.fetch_dstore("ear_height") * m;

    _msetup.message("MIND ear height:",EAR_height/cm,"cm",c);

    EAR_width = _pstore.fetch_dstore("ear_width") * m;

    _msetup.message("MIND ear width:",EAR_width/cm,"cm",c);

    
    
    VERT_x = _pstore.find_dstore("vertex_x") ?
      _pstore.fetch_dstore("vertex_x") * m : 0.0;
    VERT_y = _pstore.find_dstore("vertex_y") ?
      _pstore.fetch_dstore("vertex_y") * m : 0.0;
    VERT_z = _pstore.find_dstore("vertex_z") ?
      _pstore.fetch_dstore("vertex_z") * m : 0.0;
    

    //-------------------------------------------------------------//
    //                      | INNER DIMENSIONS |                   //
    //-------------------------------------------------------------//

    IRON_z = _pstore.fetch_dstore("passive_thickness") * cm;
    SCINT_z = _pstore.fetch_dstore("active_thickness") * cm;
    AIR_z = _pstore.fetch_dstore("air_gap") * cm;
    Al_z  = _pstore.find_dstore("bracing_thickness") ? _pstore.fetch_dstore("bracing_thickness") * cm : 0.0;
    nScint = _pstore.fetch_istore("active_layers");


    //Adjust length for integer number of pieces.
    _pieceWidth = IRON_z + nScint*(SCINT_z + Al_z) + (nScint+1)*AIR_z;
    _npieces = (int)ceil( MIND_z / _pieceWidth );
    MIND_z = _npieces * _pieceWidth;
    rel_denSI = _pstore.fetch_dstore("rel_denSI");
    rel_denAS = _pstore.fetch_dstore("rel_denAS");

    //--------------------------- VOLUMES ------------------------//
    
    MOTHER_x = MIND_x + 10 * cm; // always assume that the vertex 
    MOTHER_y = MIND_y + 10 * cm; // detector width is smaller than the MIND
    MOTHER_z = MIND_z + VERT_z + 10 * cm;  
    MOTHER_earh = EAR_height - 10 * cm;
    MOTHER_earw = EAR_width  + 10 * cm;

    // -------------------------------------------------------------//
    //                       |  MAGNETIC FIELD |                    //
    // -------------------------------------------------------------//
    
    bhep::vdouble field = _pstore.fetch_vstore("mag_field");
    BField = EVector(3,0);
    BField[0] = field[0] * tesla; BField[1] = field[1] * tesla;
    BField[2] = field[2] * tesla;

    // override constant BField if field map is present
    if(_pstore.find_sstore("mag_field_map")) {
      std::string Bmap =  _pstore.fetch_sstore("mag_field_map");
      double fieldScale = 1.0;
      if (_pstore.find_dstore("fieldScale") ) {
	fieldScale = _pstore.fetch_dstore("fieldScale");
	std::cout<<"Field Scaling is "<<fieldScale<<std::endl;
      }
      fieldScale *= IRON_z > 0 ? IRON_z/_pieceWidth : 1.0;
      BFieldMap = MINDfieldMapReader(Bmap,fieldScale);
      B_int = fieldScale * tesla;
    }
    else {
      // Default is to use a radially symmetric field 
      double fieldScale = 1.0;
      if (_pstore.find_dstore("fieldScale"))
	fieldScale = _pstore.fetch_dstore("fieldScale");
      // This uses an analytic approximation of a simulated field map 
      // std::cout<<IRON_z/_pieceWidth<<std::endl;
      fieldScale *= IRON_z > 0 ? IRON_z/_pieceWidth : 1.0;
      BFieldMap = MINDfieldMapReader(fieldScale, OctGeom, MIND_x, MIND_y);
      B_int = fieldScale * tesla;
    }
   
    //_msetup.message("Magnetic field intensity:",B_int/tesla,"tesla",c);
    
    // -------------------------------------------------------------//
    //            |  RADIATION LENGTH AND ENERGY LOSS |             //
    // -------------------------------------------------------------//
    
    X0Fe = _pstore.fetch_dstore("x0Fe") * mm;
    X0Sc = _pstore.fetch_dstore("x0Sc") * mm;
    X0AIR = _pstore.fetch_dstore("x0AIR") * m;

    double wSc = SCINT_z / (SCINT_z + AIR_z*(nScint+1)*rel_denAS);
    double X01 = (X0Sc*X0AIR) / (wSc*(X0AIR-X0Sc) + X0Sc);
    _wFe = IRON_z/(IRON_z + ((SCINT_z+AIR_z)*nScint+AIR_z)*rel_denSI*(wSc*(1-rel_denAS)+rel_denAS));
    // _wFe = IRON_z/(IRON_z + rel_denSI*(SCINT_z+rel_denAS*AIR_z));

    X0Eff = 1./(_wFe/X0Fe + wSc/X01);

    //de_dx = _pstore.fetch_dstore("de_dx") * MeV/cm;
    //de_dx = 1./(_wFe/de_dxFe + wSc/de_dxFe);

    // changed to introduce the de_dx map
    de_dx_scint = _pstore.find_dstore("de_dx_scint")?
      _pstore.fetch_dstore("de_dx_scint") * MeV/mm : 0.205 * MeV/mm;
    de_dx_min = _pstore.fetch_dstore("de_dx_min") * MeV/mm;

    _msetup.message("Radiation length:",X0Fe/cm,"cm",c);

    // -------------------------------------------------------------//
    //                       |  MEASUREMENTS |                    //
    // -------------------------------------------------------------//

    meas_dim = 2;
    meastype = _pstore.fetch_sstore("meas_type");
    
    resx = _pstore.fetch_dstore("pos_res") * cm;
    resy = _pstore.fetch_dstore("pos_res") * cm;
    
    if(_pstore.find_dstore("StepSize")){
      StepSize = _pstore.fetch_dstore("StepSize") * cm;
    }
    else{
      StepSize = 0.;
    }

}


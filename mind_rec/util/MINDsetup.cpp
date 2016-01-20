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

  delete mother;

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
  //Volume* mother;
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

  //dict::Key vol_name = "Detector";
  //const dict::Key vol_name = "Detector";
  const dict::Key vert_name = "VertDetector";
  // Volume* det = new Box(pos,xaxis,yaxis,MIND_x/2,MIND_y/2,MIND_z/2);
    
  dict::Key vol_name;

  for(map<string, std::vector<double> >::const_iterator it = _gdml_pos_map.begin();
      it != _gdml_pos_map.end(); ++it)
    {
      std::map<string,std::vector<double> >::iterator it_int;
      vol_name = it->first;
      pos[0]= it->second[0];
      pos[1]= it->second[1];
      pos[2]= it->second[2];

      if(isdigit(it->first.at(it->first.length() -1)))
	{
	  //If the string ends on a digit, remove it before using as a key
	  it_int =_gdml_solid_map.find(it->first.substr(0,it->first.length()-1));
	} 
      else
	{
	  it_int =_gdml_solid_map.find(it->first);
	}
      
      if(it_int!= _gdml_solid_map.end())
	{
	  //cout << it_int->first<<endl;
	  //det = new Box(pos, zaxis, xaxis, it_int->second[2]/2.,
	  //		it_int->second[0]/2, it_int->second[1]/2);

	  detVector.push_back(new Box(pos, zaxis, xaxis, it_int->second[2]/2.,
	  		      it_int->second[0]/2, it_int->second[1]/2));

	  //_gsetup.add_volume("mother",vol_name,det);
	  _gsetup.add_volume("mother",vol_name,detVector.back());

	  vector<double> tempVector;
	  tempVector.push_back(pos[2]);
	  tempVector.push_back(it_int->second[2]/2);

	  _moduleDataMap[vol_name] = tempVector;
	}
      //  std::cout << it->first << " " << it->second[0] << " " << it->second[1] << " "  <<it->second[2]<< "\n";
    }
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
  
  Fe_slab = new MINDplate(fe_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,
				  IRON_z/2, EAR_width, EAR_height);
  
  _gsetup.add_volume(det_vol,iron_name,Fe_slab);
  
  _gsetup.set_volume_property(iron_name,"X0",X0Fe);
  
  /*
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
  */
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

  _zaxis = EVector(3,0);
  _zaxis[2]=1;

  string current_string;
  dict::Key vol_name;

  _generalBFieldMap = MINDfieldMapReader(Bmap,_fieldScale);
 
  /*

  Fill the detector subvolumes with the correct properties by iterating over each part, finding
  the ammount of iron and scintilator and the total size of the module.
  */
  for(map<string, std::vector<double> >::const_iterator it = _gdml_pos_map.begin();
      it != _gdml_pos_map.end(); ++it)
    {

      double fieldScale = _fieldScale;
      //double fieldScale = 1;
      string name = it->first;
      int numScint=0;
      int numFe=0;
      vol_name = it->first;
      
      // From the name, count the number of S and F.
      for(unsigned int len= 0; len < name.length(); len++)
	{
	  current_string = name.at(len);
	  if("S" == current_string)
	    {
	      numScint++;
	    }
	  else if("F" == current_string)
	    {
	      numFe++;
	    }
	  else if("M" == current_string)
	    {
	      if(numFe ==3)
		{
		  numFe *= 6;
		  numScint *= 6;
		}
	      else
		{
		  numFe *=3;
		  numScint *=3;
		}
	      break;
	    }
	  else
	    {
	      break;
	    }
	}

      // Find the solid reference for the module
      std::map<string,std::vector<double> >::iterator it_int;
      
      if(isdigit(it->first.at(it->first.length() -1)))
	{
	  //If the string ends on a digit, remove it before using as a key
	  it_int =_gdml_solid_map.find(it->first.substr(0,it->first.length()-1));
	} 
      else
	{
	  it_int =_gdml_solid_map.find(it->first);
	}
  
	if ((numScint != 0 || numFe != 0) & it_int!= _gdml_solid_map.end())
	{
	  // Calculate the appropriate properties and add them to the boxes
	  //cout<<"YES"<<endl;
	  
	  // Should be done in readParam
	  double wSc = SCINT_z / (SCINT_z + AIR_z*(numScint+1)*rel_denAS);
	  double X01 = (X0Sc*X0AIR) / (wSc*(X0AIR-X0Sc) + X0Sc);
	  double wFe = IRON_z/(IRON_z + ((SCINT_z+AIR_z)*numScint+AIR_z)*rel_denSI*(wSc*(1-rel_denAS)+rel_denAS));
	  
	  _wFe = wFe;

	  double X0Eff = 1./(wFe/X0Fe + wSc/X01);
	  //double X0Eff = 1./(wFe/X0Fe + wSc/X0Sc);
	  X0EffVec.push_back(X0Eff);

	  //double length = numScint * AIR_z +numScint * SCINT_z + numFe * IRON_z;
	  double length = it_int->second[2]; // Simply taken from the solid reference.

	  double de_dx = (numScint * SCINT_z * de_dx_scint + numFe * IRON_z * de_dx_fe)/length;

	  _moduleDataMap[vol_name].push_back(de_dx);

	  std::cout<<"modulelength "<<length<<" numFe "<<numFe<<" numScint "<<numScint<<std::endl;
	  std::cout<<"iron_z "<<IRON_z<<" scint_z "<<SCINT_z<<std::endl;

	  fieldScale *= numFe > 0 ? IRON_z*numFe/length : 0;
	  //fieldScale *= IRON_z *numFe;
	  
	  std::cout<<"Local Field Scaling is "<<fieldScale<<std::endl;
	  //BFieldMap = MINDfieldMapReader(Bmap,fieldScale);
	  BFieldMapVec.push_back(new MINDfieldMapReader(Bmap,fieldScale));
	  //BFieldMap = MINDfieldMapReader(Bmap,fieldScale);
	  // Let the scale always be one and fill the proper one in moduleDataMap, the scaling is done in 
	  // getBfield below.

	  double copy = fieldScale;
	  _moduleDataMap[vol_name].push_back(copy);

	  //_de_dx_map = new DeDxMap(de_dx*MeV/mm);
	  de_dx_map_vec.push_back(new DeDxMap(de_dx*MeV/mm));
	  //_de_dx_map_scint = new DeDxMap(de_dx_scint*MeV/mm);

	  //_gsetup.set_volume_property(vol_name,RP::de_dx_map,*_de_dx_map);
	  _gsetup.set_volume_property(vol_name,RP::de_dx_map,*de_dx_map_vec.back());
	  //_gsetup.set_volume_property(vol_name,"X0",X0Eff);
	  _gsetup.set_volume_property(vol_name,"X0",X0EffVec.back());
	  //_gsetup.set_volume_property_to_sons("mother",RP::BFieldMap,BFieldMap);
	  //_gsetup.set_volume_property(vol_name,RP::BFieldMap,BFieldMap);
	  _gsetup.set_volume_property(vol_name,RP::BFieldMap,*BFieldMapVec.back());

	}
    }  
}


EVector MINDsetup::getBField(EVector pos){
  // Has become more advance now that we have subdetectors, find the fieldScale in moduelDataMap 
  //and return the value properly scaled.

  // Return the value from a general fieldmap scaled correctly for the subdetector.

  EVector BfieldVector = _generalBFieldMap.vector(pos);
  
  double properScale = 0;
  double zCoord = pos[2];

  for (std::map<dict::Key,vector<double> >::iterator it=_moduleDataMap.begin();
       it!=_moduleDataMap.end(); ++it)
    {
      // vol_name, module position z, module size z, wFe, magfieldScale
      double module_pos = it->second[0];
      double module_half_size = it->second[1];
      double fieldScale = it->second[3];
      
      if(zCoord>=(module_pos - module_half_size) && zCoord<=(module_pos + module_half_size))
	{
	  properScale = fieldScale;
	  break;
	}

    }



  return properScale*BfieldVector;
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
    
    _fieldScale = 1.0;
    //if (_pstore.find_dstore("fieldScale") ) {
    //_fieldScale = _pstore.fetch_dstore("fieldScale");
      //std::cout<<"Field Scaling is "<<_fieldScale<<std::endl;
    //}

    Bmap = _pstore.fetch_sstore("mag_field_map");
	  

    /*
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
    */
    //_msetup.message("Magnetic field intensity:",B_int/tesla,"tesla",c);
    
    // -------------------------------------------------------------//
    //            |  RADIATION LENGTH AND ENERGY LOSS |             //
    // -------------------------------------------------------------//
    
    X0Fe = _pstore.fetch_dstore("x0Fe") * mm;
    X0Sc = _pstore.fetch_dstore("x0Sc") * mm;
    X0AIR = _pstore.fetch_dstore("x0AIR") * mm;

    de_dx_scint = _pstore.fetch_dstore("de_dx_scint") * MeV/mm;
    de_dx_fe = _pstore.fetch_dstore("de_dx_fe") * MeV/mm;

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


    // -------------------------------------------------------------//
    //                  |  USE THE GDML PARSED FILE    |            //
    // -------------------------------------------------------------//

  //read gdml_parsed file
  _gdml_parsed_path = _pstore.fetch_sstore("xml_parsed");


  // Fill _gdml_solid_map with all the solid references and
  // Fill  _gdml_pos_map with all the actual parts with their positions.

  std::ifstream file;
  file.open (_gdml_parsed_path.c_str());

  string word;
  string temp;
  string type;
  double x;
  double y;
  double z;
  std::vector<double> temp_vector;
  while(file >> type)
    {
      if (type == "SOLID")
	{
	  file >> word;
	  temp_vector.clear();
	  file >> temp;
	  x = mm*atof(temp.c_str());
	  temp_vector.push_back(x);
	  file >> temp;
	  y = mm*atof(temp.c_str());
	  temp_vector.push_back(y);
	  file >> temp;
	  z = mm*atof(temp.c_str());
	  temp_vector.push_back(z);
	  _gdml_solid_map[word] = temp_vector;
	}
      else if (type == "POS")
	{
	  file >> word;
	  temp_vector.clear();
	  file >> temp;
	  x = mm*(double)atof(temp.c_str());
	  temp_vector.push_back(x);
	  file >> temp;
	  y = mm*(double)atof(temp.c_str());
	  temp_vector.push_back(y);
	  file >> temp;
	  z = mm*(double)atof(temp.c_str());
	  temp_vector.push_back(z);
	  _gdml_pos_map[word] = temp_vector;
	}
    

    }


}


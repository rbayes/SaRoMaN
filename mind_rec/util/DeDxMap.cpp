#include "DeDxMap.h"
#include <recpack/Definitions.h>

//*******************************************************
Recpack::DeDxMap::DeDxMap(double dedx_min):DynamicProperty(RP::de_dx_map) {
//*******************************************************

  _de_dx_min = dedx_min;

   //  _de_dx_min_def = 0.195*MeV/unit::mm;

  _melectron = 0.51099906 * MeV;
  _mproton = 938.27231 * MeV;
  _mmuon = 105.65837 * MeV;
  _mpion = 139.57018 * MeV;
  _mkaon = 493.677 * MeV;

  // map for muons in scintillator

  _pbins[0]=0.;               _de_dx[0] =18.44*MeV/mm;
  _pbins[1]=10*MeV;     _de_dx[1] =9.19*MeV/mm;  
  _pbins[2]=20*MeV;     _de_dx[2] =3.16*MeV/mm;
  _pbins[3]=30*MeV;     _de_dx[3] =1.72*MeV/mm;
  _pbins[4]=40*MeV;     _de_dx[4] =1.10*MeV/mm;
  _pbins[5]=50*MeV;     _de_dx[5] =0.78*MeV/mm;
  _pbins[6]=60*MeV;     _de_dx[6] =0.60*MeV/mm;
  _pbins[7]=70*MeV;     _de_dx[7] =0.49*MeV/mm;
  _pbins[8]=80*MeV;     _de_dx[8] =0.42*MeV/mm;
  _pbins[9]=90*MeV;     _de_dx[9] =0.36*MeV/mm;
  _pbins[10]=100*MeV;   _de_dx[10]=0.31*MeV/mm;
  _pbins[11]=120*MeV;   _de_dx[11]=0.27*MeV/mm;
  _pbins[12]=140*MeV;   _de_dx[12]=0.244*MeV/mm;
  _pbins[13]=160*MeV;   _de_dx[13]=0.228*MeV/mm;
  _pbins[14]=180*MeV;   _de_dx[14]=0.217*MeV/mm;
  _pbins[15]=200*MeV;   _de_dx[15]=0.202*MeV/mm;
  _pbins[16]=300*MeV;   _de_dx[16]=0.195*MeV/mm;
  _pbins[17]=400*MeV;   _de_dx[17]=0.196*MeV/mm;
  _pbins[18]=500*MeV;   _de_dx[18]=0.199*MeV/mm;
  _pbins[19]=600*MeV;   _de_dx[19]=0.202*MeV/mm;
  _pbins[20]=700*MeV;   _de_dx[20]=0.206*MeV/mm;
  _pbins[21]=800*MeV;   _de_dx[21]=0.209*MeV/mm;
  _pbins[22]=900*MeV;   _de_dx[22]=0.212*MeV/mm;
  _pbins[23]=1000*MeV;  _de_dx[23]=0.224*MeV/mm;
  _pbins[24]=2000*MeV;  _de_dx[24]=0.242*MeV/mm;
  _pbins[25]=3000*MeV;  _de_dx[25]=0.253*MeV/mm;
  _pbins[26]=4000*MeV;  _de_dx[26]=0.261*MeV/mm;
  _pbins[27]=5000*MeV;  _de_dx[27]=0.274*MeV/mm;


  // map for electrons in scintillator

  _de_dx_elec[0]=0.240*MeV/mm;
  _de_dx_elec[1]=0.249*MeV/mm;  
  _de_dx_elec[2]=0.266*MeV/mm;
  _de_dx_elec[3]=0.276*MeV/mm;
  _de_dx_elec[4]=0.283*MeV/mm;
  _de_dx_elec[5]=0.289*MeV/mm;
  _de_dx_elec[6]=0.293*MeV/mm;
  _de_dx_elec[7]=0.297*MeV/mm;
  _de_dx_elec[8]=0.299*MeV/mm;
  _de_dx_elec[9]=0.302*MeV/mm;

  _de_dx_elec[10]=0.304*MeV/mm;
  _de_dx_elec[11]=0.307*MeV/mm;
  _de_dx_elec[12]=0.309*MeV/mm;
  _de_dx_elec[13]=0.311*MeV/mm;
  _de_dx_elec[14]=0.312*MeV/mm;


  _de_dx_elec[15]=0.315*MeV/mm;
  _de_dx_elec[16]=0.318*MeV/mm;
  _de_dx_elec[17]=0.319*MeV/mm;
  _de_dx_elec[18]=0.320*MeV/mm;
  _de_dx_elec[19]=0.321*MeV/mm;
  _de_dx_elec[20]=0.321*MeV/mm;
  _de_dx_elec[21]=0.321*MeV/mm;
  _de_dx_elec[22]=0.322*MeV/mm;

  _de_dx_elec[23]=0.322*MeV/mm;
  _de_dx_elec[24]=0.323*MeV/mm;
  _de_dx_elec[25]=0.323*MeV/mm;
  _de_dx_elec[26]=0.323*MeV/mm;
  _de_dx_elec[27]=0.323*MeV/mm;

}


//*******************************************************
double Recpack::DeDxMap::property(const State& state) const{
//*******************************************************


  double de_dx = 0;

  double qop = state.vector()[state.dim()-1];
  
  double p=10000;  // default momentum
  if (qop !=0) 
    p = 1/(fabs(qop));
  
  Key PID= "";
  if (state.names().has_key(RP::PID)) 
    PID = state.name(RP::PID);

  // Use momentum based on particle type

  double p_pid = p;
  /*
  if (PID==RP::muon)
    p_pid = p;
  else if (PID==RP::electron)
    p_pid = p;
  else if (PID==RP::pion)
    p_pid = p*_mpion/_mmuon;
  else if (PID==RP::proton)
    p_pid = p*_mproton/_mmuon;
  else if (PID==RP::kaon)
    p_pid = p*_mkaon/_mmuon;
  */

  // if none of this Hypothesis it will assume the Muon
  if (PID=="Muon")
    p_pid = p;
  else if (PID=="Electron")
    p_pid = p;
  else if (PID=="Pion")
    p_pid = p*_mpion/_mmuon;
  else if (PID=="Proton")
    p_pid = p*_mproton/_mmuon;
  else if (PID=="Kaon")
    p_pid = p*_mkaon/_mmuon;

  
  // Find the bin number
  int pbin=0;
  if (p_pid>_pbins[NPBINS-1]) pbin = NPBINS-1;
  else{
    for (int i=0;i<NPBINS-1;i++){
      if (p_pid>_pbins[i] && p_pid<_pbins[i+1]){
	pbin = i;
	break;
      }
    }
  }
  



  if (PID==RP::electron)
    de_dx = - _de_dx_elec[pbin]*_de_dx_min/_de_dx_min_def;
  else
    de_dx = - _de_dx[pbin]*_de_dx_min/_de_dx_min_def;

  if (verbosity(VERBOSE))
    std::cout << "DeDxMap: pid=" << PID << ", p=" << p << ", p'=" << p_pid << ", pbin=" << pbin << ", p0=" << _pbins[pbin] << " --> " << de_dx << std::endl;

  return de_dx;

}
///
/*/*************************************************************
double Recpack:: DeDxMap :: getDeDx_map(double mom){
  //*************************************************************
  double mapped_de_dx =0;

  if (mom > _pbins[NPBINS]){
    mapped_de_dx = _de_dx[NPBINS];
   
  }
 
  for (int i=0; i<=NPBINS; i++){
    if (mom > _pbins[i] && mom < _pbins[i+1])
      mapped_de_dx = _de_dx[i];
    break;
  }
      
  return mapped_de_dx;
  }*/

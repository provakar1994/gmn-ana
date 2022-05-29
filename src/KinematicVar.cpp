#include "../include/KinematicVar.h"

namespace kine {

  double Q2(double ebeam, double eeprime, double etheta){
    return 2.0*ebeam*eeprime*(1.0-cos(etheta));
  }
  //--------------------------------------------
  double W2(double ebeam, double eeprime, double Q2){
    return pow(constant::Mp,2.0) + 2.0*constant::Mp*(ebeam-eeprime) - Q2;
  }
  //--------------------------------------------
  double W(double ebeam, double eeprime, double Q2){
    if(kine::W2(ebeam,eeprime,Q2)>0)
      return sqrt(kine::W2(ebeam,eeprime,Q2));
    else
      return 0;
  }
  //--------------------------------------------

} //::kine

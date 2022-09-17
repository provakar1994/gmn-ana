#ifndef EXP_CONSTANTS_H
#define EXP_CONSTANTS_H

#include "TMath.h"
// #include "SBSconfig.h"

namespace expconst {

  // Detector variables
  static const int shcol = 7;
  static const int shrow = 27;
  static const int pscol = 2;
  static const int psrow = 26;
  static const int hcalcol = 12;
  static const int hcalrow = 24;
  
  // Constant for the entire experiment
  // target
  static const double tarlen = 15.0;  //cm
 
  // LH2
  static const double lh2tarrho = 0.0723;     //g/cc, target density
  static const double lh2cthick = 0.02;       //cm, target cell thickness
  static const double lh2uwallthick = 0.0145; //cm, upstream wall thickness
  static const double lh2dwallthick = 0.015;  //cm, downstream wall thickness

  // LD2
  static const double ld2tarrho = 0.169;      //g/cc, target density

  // magnet
  static const double dipolegap = 1.22;  //m

  // shieldling
  static const double Alrho = 2.7; //g/cc

  // Following quantities vary with configuration
  double ebeam(int config);     //GeV
  double bbtheta(int config);   //deg
  double bbdist(int config);    //m
  double sbstheta(int config);  //deg
  double sbsdist(int config);   //m
  double hcaldist(int config);  //m
}

#endif

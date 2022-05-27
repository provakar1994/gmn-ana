#ifndef EXP_CONSTANTS_H
#define EXP_CONSTANTS_H

namespace expconst {

  // Detector variables
  static const int shcol = 7;
  static const int shrow = 27;
  static const int pscol = 2;
  static const int psrow = 26;
  static const int hcalcol = 12;
  static const int hcalrow = 24;
  
  // Constant for the entire experiment
  static const double dipolegap = 1.22;  //m

  // Following quantities vary with configuration
  double ebeam(int config);  //GeV
  double bbtheta(int config);  //deg
  double bbdist(int config);  //m
  double sbstheta(int config);  //deg
  double sbsdist(int config);  //m
  double hcaldist(int config);  //m
  
}

#endif

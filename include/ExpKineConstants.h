#ifndef EXP_CONSTANTS_H
#define EXP_CONSTANTS_H

namespace expconst {

  // Constant for the entire experiment
  double dipolegap = 1.22;  //m

  // Following quantities vary with configuration
  double ebeam(int config);  //GeV
  double bbtheta(int config);  //deg
  double bbdist(int config);  //m
  double sbstheta(int config);  //deg
  double sbsdist(int config);  //m
  double hcaldist(int config);  //m
  
}

#endif

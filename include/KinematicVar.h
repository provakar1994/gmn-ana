#ifndef KINE_VAR_H
#define KINE_VAR_H

#include <cmath>
#include "../include/Constants.h"

namespace kine{

  double Q2(double ebeam, double eeprime, double etheta); // GeV, GeV, rad
  double W2(double ebeam, double eeprime, double Q2); // GeV, GeV, GeV2
  double W(double ebeam, double eeprime, double Q2); // GeV, GeV, GeV2
  
}

#endif

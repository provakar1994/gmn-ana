#ifndef UTIL_CALIB_COEFF_H
#define UTIL_CALIB_COEFF_H

// a struct to store calibration coefficients 
#include <cstdlib>
#include <string> 

typedef struct calibCoeff {
   std::string dev;      // name of device: unser, u1, unew, d1, d3, d10, dnew
   std::string info;     // user info 
   double pedestal;      // pedestal (beam off) [Hz] 
   double pedestalErr;   // pedestal (beam off) error [Hz]
   double offset;        // fitted offset [Hz]  
   double offsetErr;     // fitted offset error [Hz]   
   double slope;         // fitted slope (gain)  [Hz/uA] 
   double slopeErr;      // fitted slope (gain) error [Hz/uA] 
   int runMin;           // run range: minimum    
   int runMax;           // run range: maximum   
   calibCoeff(): 
      dev("NONE"),info("NONE"),
      pedestal(0),pedestalErr(0),
      offset(0),offsetErr(0),
      slope(0),slopeErr(0),
      runMin(-1),runMax(-1)
   {}
   // useful functions 
   void Print(){
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "dev:        " << dev                 << std::endl;
      std::cout << "run range:  " << runMin << ", " << runMax << std::endl;
      std::cout << "pedestal:   " << Form("%.3lf ± %.3lf",pedestal,pedestalErr)  << std::endl;
      std::cout << "offset:     " << Form("%.3lf ± %.3lf",offset,offsetErr)      << std::endl;
      std::cout << "slope:      " << Form("%.3lf ± %.3lf",slope,slopeErr)        << std::endl;
      std::cout << "----------------------------------------" << std::endl;
   }
} calibCoeff_t; 

#endif 

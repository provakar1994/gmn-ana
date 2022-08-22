#ifndef UTIL_PRODUCED_VARIABLE_H 
#define UTIL_PRODUCED_VARIABLE_H 

// a simple struct to hold produced data  
#include <cstdlib>
#include <iostream>
#include <string>
#include "TString.h"

typedef struct producedVariable {
   std::string dev;            // device name 
   std::string beam_state;     // beam state (on or off)
   double time;                // timestamp  
   double mean;                // mean 
   double stdev;               // stdev 
   int group;                  // group/categorization 
   // constructor 
   producedVariable():
      dev("NONE"),beam_state("NONE"),time(0),mean(0),stdev(0),group(0) 
   {}
   // useful functions 
   void Print(){
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "dev:        " << dev                 << std::endl;
      std::cout << "beam_state: " << beam_state          << std::endl;
      std::cout << "group:      " << group               << std::endl;
      std::cout << "time:       " << Form("%.3lf",time)  << std::endl;
      std::cout << "mean:       " << Form("%.3lf",mean)  << std::endl;
      std::cout << "stdev:      " << Form("%.3lf",stdev) << std::endl;
      std::cout << "----------------------------------------" << std::endl;
   } 
} producedVariable_t; 

#endif 

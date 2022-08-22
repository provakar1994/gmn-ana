#ifndef UTIL_CHARGE_DF
#define UTIL_CHARGE_DF

// a struct to collect charge calculation results

#include <cstdlib>
#include <iostream> 

typedef struct charge { 
   std::string info;   // general run information
   double totalTime;   // total time of run (sec)  
   double value;       // charge in Coulombs 
   double error;       // uncertainty in Coulombs
   int runNumber;      // run number 
   // constructor 
   charge(): 
      info("NONE"),totalTime(0),value(0),error(0),runNumber(0)
   {} 
   // useful functions 
   void Print(){
      std::cout << "-------------------------" << std::endl;
      std::cout << "Run number: " << runNumber << std::endl; 
      std::cout << "Total time: " << Form("%.3lf sec",totalTime)     << std::endl;
      std::cout << "Charge:     " << Form("%.3E ± %.3E",value,error) << std::endl;
      std::cout << "-------------------------" << std::endl;
   }
} charge_t; 

#endif 

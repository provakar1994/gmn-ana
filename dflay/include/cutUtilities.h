#ifndef CUT_UTIL_H
#define CUT_UTIL_H 

// a namespace for helpful functions useful in 
// BCM data analysis 

#include <cstdlib>
#include <iostream>
#include <vector> 

#include "cut.h"
#include "scalerData.h"
#include "epicsData.h"

namespace cut_util {
   int LoadCuts(const char *inpath,std::vector<cut_t> &data);
   int LoadCuts_json(const char *inpath,std::vector<cut_t> &data); 
   int LoadCuts_epics_json(const char *inpath,std::vector<cut_t> &data);
   int ApplyCut(std::string var,double lo,double hi,std::vector<scalerData_t> in,std::vector<scalerData_t> &out); 
   int ApplyCut(std::string var,double lo,double hi,std::vector<epicsData_t> in,std::vector<epicsData_t> &out); 
   int ApplyCuts(double cutLo,double cutHi,std::vector<double> x,std::vector<double> y,std::vector<double> &out); 
   int ApplyCuts(std::vector<cut_t> cutList,std::vector<scalerData_t> in,std::vector<scalerData_t> &out); 
   int ApplyCuts_alt(std::vector<cut_t> cutList,std::vector<scalerData_t> in,std::vector<scalerData_t> &out); 
   int ApplyCuts_alt(std::vector<cut_t> cutList,std::vector<epicsData_t> in,std::vector<epicsData_t> &out); 
   int GetStatsWithCuts(std::vector<double> x,std::vector<double> y,
                        double cutLo,double cutHi,double &mean,double &stdev);

   int GetStatsWithCuts(double cutLo,double cutHi,
                        std::vector<double> x,std::vector<double> y1,std::vector<double> y2,
                        double &mean1,double &stdev1,double &mean2,double &stdev2);  

   int ApplyCut(cut_t aCut,std::vector<scalerData_t> data,std::vector<scalerData_t> &out);  

   bool PassCut(cut_t cut,scalerData_t event); 
   bool PassCut_alt(cut_t cut,scalerData_t event); 
   bool PassCut_alt(cut_t cut,epicsData_t event); 
}

#endif 

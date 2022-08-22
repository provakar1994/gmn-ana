#ifndef UTIL_BCM_H
#define UTIL_BCM_H 

// a namespace for helpful functions useful in 
// beam charge/current monitor (BCM) data analysis  

#include <cstdlib>
#include <iostream>
#include <vector> 

#include "cut.h"
#include "charge.h"
#include "codaRun.h"
#include "calibCoeff.h"
#include "scalerData.h"
#include "beamData.h"
#include "producedVariable.h"

namespace bcm_util {
   int Print(producedVariable_t data); 
   int GetData(std::string var,std::vector<producedVariable_t> data,std::vector<double> &x,std::vector<double> &dx); 
   int WriteToFile(const char *outpath,std::vector<producedVariable_t> data); 
   int WriteToFile_cc(const char *outpath,std::vector<calibCoeff_t> data);
   int LoadCalibrationCoefficients(const char *inpath,std::vector<calibCoeff_t> &data); 
   int LoadFittedOffsetGainData(const char *inpath,std::vector<calibCoeff_t> &data); 
   int LoadPedestalData(const char *inpath,std::vector<calibCoeff_t> &data); 
   int LoadProducedVariables(const char *inpath,std::vector<producedVariable_t> &data);
   int LoadConfigPaths(const char *inpath,std::vector<std::string> &label,std::vector<std::string> &path); 
   int LoadRuns(const char *inpath,std::vector<codaRun_t> &runList); 
   int LoadCuts(const char *inpath,std::vector<cut_t> &data);
   int ApplyCuts(double cutLo,double cutHi,std::vector<double> x,std::vector<double> y,std::vector<double> &out); 
   int GetStatsWithCuts(std::vector<double> x,std::vector<double> y,
                        double cutLo,double cutHi,double &mean,double &stdev);
   int ConvertToCurrent(calibCoeff_t cc,std::vector<producedVariable_t> unser_ps,
	                std::vector<producedVariable_t> &unser_cur);

   int CalculateStatsForBeamState(std::string beamState,std::vector<producedVariable_t> data,std::vector<producedVariable_t> &out,
                                  double &MEAN,double &STDEV,std::string LOG_PATH); 
   int CalculatePedestalSubtraction(std::vector<producedVariable_t> data,std::vector<producedVariable_t> &out,
                                    std::string LOG_PATH="NONE",std::string PLOT_PATH="NONE"); 
   int SubtractBaseline(std::vector<producedVariable_t> on,std::vector<producedVariable_t> off, 
                        std::vector<producedVariable_t> &diff,bool isDebug=false);

   int GetCharge(std::string var,std::vector<scalerData_t> runData,charge_t &out);
   double GetChargeError(double Q,double I, double dI,double t,double dt); 

   int GetStats_byRun(std::string var,std::vector<scalerData_t> data,
                      std::vector<double> &RUN,std::vector<double> &MEAN,std::vector<double> &STDEV); 
   int GetStats_byRun_epics(std::string var,std::vector<epicsData_t> data,
                            std::vector<double> &RUN,std::vector<double> &MEAN,std::vector<double> &STDEV); 

   TGraph *GetTGraph(std::string xAxis,std::string yAxis,std::vector<epicsData_t> data); 
   TGraph *GetTGraph(std::string xAxis,std::string yAxis,std::vector<scalerData_t> data); 
   TGraph *GetTGraph_singleRun(int run,std::string xAxis,std::string yAxis,std::vector<scalerData_t> data); 
   TGraphErrors *GetTGraphErrors_byRun(std::string var,std::vector<scalerData_t> data);
   TGraphErrors *GetTGraphErrors_byRunByUnserCurrent(std::string var,std::vector<scalerData_t> data); 
}

#endif 

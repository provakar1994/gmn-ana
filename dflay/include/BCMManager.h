#ifndef BCM_MANAGER_H
#define BCM_MANAGER_H

// a class that can produce plots of the BCM variables 
// as a function of event number or time in the run 

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>  

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TDatime.h"

#include "THaRun.h"

#include "scalerData.h"
#include "epicsData.h"
#include "calibCoeff.h"

class BCMManager {

   private:
      bool fIsDebug,fCalculateCurrent;
      int fVerbosity;  
      int fLastRun,fLastRunEvtCntLeft,fLastRunEvtCntSBS,fLastRunEvtCntEPICS;
      int fEvtCntrLeft,fEvtCntrSBS,fEvtCntrEPICS;
      double fLastTimeLeft,fLastTimeSBS;

      unsigned long fUTCTimeStamp; 
      std::string fTimeStamp; 

      // calibration coefficients to convert rate to current [uA]
      // old style: unser, upstream and downstream BCMs in one vector 
      std::vector<calibCoeff_t> fCC;   

      // calibration coefficients for multiple run ranges 
      std::vector<calibCoeff_t> fccUnser;
      std::vector<calibCoeff_t> fccU1,fccUnew;
      std::vector<calibCoeff_t> fccD1,fccD3,fccD10,fccDnew;

      // data structs to hold scaler data  
      std::vector<scalerData_t> fLeft,fSBS;
      std::vector<epicsData_t> fEPICS;

      std::vector<int> fRunList; // list of runs read into the manager  

      bool IsBad(double v);  
      int CheckFile(const char *filePath); 
      int LoadCalibrationCoefficients(const char *type,const char *filePath);
      int GetCalibrationCoeff(int run,std::string dev,calibCoeff_t &data); 
      int ApplyCalibrationCoeff(scalerData_t &data); 

   public: 
      BCMManager(const char *filePath="NONE",const char *ccDirPath="NONE",bool isDebug=false);
      ~BCMManager(); 

      void SetVerbosity(int v)           { fVerbosity        = v; } 
      void SetDebug(bool v=true)         { fIsDebug          = v; }
      void CalculateCurrent(bool v=true) { fCalculateCurrent = v; }  
      void Print(const char *arm);
      void Clear(); 

      int GetTimeStamp(TFile *f,std::string &timeStamp,unsigned long &utc); 
      int LoadFile(const char *filePath,int runNumber=0);
      int LoadDataFromTree(const char *filePath,const char *treeName,int runNumber=0);
      int LoadEPICSDataFromTree(const char *filePath,int runNumber=0);
      int LoadCalibrationCoefficients(const char *dirName); 

      int GetVector(const char *arm,const char *var,std::vector<double> &v); 
      int GetVector_scaler(const char *arm,std::vector<scalerData_t> &data); 
      int GetVector_scaler(const char *arm,int run,std::vector<scalerData_t> &data); 
      int GetVector_epics(std::vector<epicsData_t> &data);
      int GetVector_epics(int run,std::vector<epicsData_t> &data);

      int GetRunList(std::vector<int> &r){
	 const int NR = fRunList.size();
	 for(int i=0;i<NR;i++) r.push_back(fRunList[i]);
	 return 0; 
      } 

      TH1F * GetTH1F(const char *arm,const char *var_name,int NBin,double min,double max);

      TH2F * GetTH2F(const char *arm,const char *xAxis,const char *yAxis,
	             int NXBin,double xMin,double xMax,int NYBin,double yMin,double yMax);  

      TGraph * GetTGraph(const char *arm,const char *xAxis,const char *yAxis);  

}; 

#endif 

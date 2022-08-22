#ifndef BEAM_MANAGER_H
#define BEAM_MANAGER_H

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

#include "beamData.h"
#include "epicsData.h"
#include "calibCoeff.h"

class BeamManager {

   private:
      bool fIsDebug,fCalculateCurrent;
      int fVerbosity;  
      int fEvtCntrLeft,fEvtCntrSBS,fEvtCntrEPICS;
      double fLastTimeLeft,fLastTimeSBS;

      double fS_bpmA;      // position of BPM A relative to the pivot 
      double fDeltaS_bpm;  // distance between BPMA and BPMB  

      std::vector<int> fRunList; 

      // calibration coefficients for multiple run ranges 
      // std::vector<calibCoeff_t> fccUnser;
      // std::vector<calibCoeff_t> fccU1,fccUnew;
      // std::vector<calibCoeff_t> fccD1,fccD3,fccD10,fccDnew;

      // data structs to hold scaler data  
      std::vector<beamData_t> fLeft,fSBS;
      std::vector<epicsData_t> fEPICS;

      bool IsBad(double v);  
      int CheckFile(const char *filePath); 
      // int LoadCalibrationCoefficients(const char *type,const char *filePath);
      // int GetCalibrationCoeff(int run,std::string dev,calibCoeff_t &data); 
      // int ApplyCalibrationCoeff(beamData_t &data); 

   public: 
      BeamManager(const char *filePath="NONE",const char *ccFilePath="NONE",bool isDebug=false);
      ~BeamManager(); 

      void SetVerbosity(int v)           { fVerbosity        = v; } 
      void SetDebug(bool v=true)         { fIsDebug          = v; }
      void Print(const char *arm);
      void Clear(); 

      int LoadFile(const char *filePath,int runNumber=0);
      int LoadDataFromTree(const char *filePath,const char *treeName,int runNumber=0);
      int LoadEPICSDataFromTree(const char *filePath,int runNumber=0);
      // int LoadCalibrationCoefficients(const char *dirName); 

      int GetVector(const char *arm,const char *var,std::vector<double> &v); 
      int GetVector_beam(const char *arm,std::vector<beamData_t> &data); 
      int GetVector_epics(std::vector<epicsData_t> &data);

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

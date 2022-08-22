// apply cuts to BCM data
// - input: run list and cut list
// - output: mean and stdev of results of defined cuts   

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"

#include "./include/bcm.h"
#include "./include/codaRun.h"
#include "./include/producedVariable.h"
#include "./src/CSVManager.cxx"
#include "./src/JSONManager.cxx"
#include "./src/ABA.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

int GetStats(std::string var,std::vector<scalerData_t> data,std::vector<double> &out); 
int GetStats(std::string var,std::vector<epicsData_t> data,std::vector<double> &out); 

int bcmCalibProcessCuts(const char *confPath,const char *bcmName){
 
   int rc=0;
   std::string devName = bcmName; 

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix  = jmgr->GetValueFromKey_str("ROOTfile-path");
   std::string tag     = jmgr->GetValueFromKey_str("tag");
   delete jmgr; 

   char run_path[200],cut_path[200],cut_path_e[200],out_dir[200],log_path[200];
   sprintf(run_path  ,"./input/%s/runlist.csv"           ,tag.c_str());
   sprintf(cut_path  ,"./input/%s/cuts/cutlist.csv"      ,tag.c_str());
   sprintf(cut_path_e,"./input/%s/cuts/cutlist_epics.csv",tag.c_str());
   sprintf(out_dir   ,"./output/%s"                      ,tag.c_str()); 
   sprintf(log_path  ,"./output/%s/log/process-cuts.txt" ,tag.c_str()); 

   std::string msg = "Processing cuts for variable: " + devName; 

   util_df::LogMessage(log_path,"========== BCM CALIBRATION PROCESS CUTS ==========",'a'); 
   util_df::LogMessage(log_path,msg.c_str(),'a'); 

   // make output directory (done by python code!) 
   // util_df::MakeDirectory(out_dir);  

   // load BCM data
   BCMManager *mgr = new BCMManager("NONE","NONE",false);

   std::vector<codaRun_t> runList;
   rc = util_df::LoadRunList(run_path,prefix.c_str(),runList);
   if(rc!=0) return 1;

   util_df::LoadBCMData(runList,mgr);
   if(rc!=0) return 1;

   // load cuts 
   std::vector<cut_t> cutList; 
   rc = bcm_util::LoadCuts(cut_path,cutList); 
   if(rc!=0) return 1; 

   std::vector<cut_t> cutListEPICS; 
   rc = bcm_util::LoadCuts(cut_path_e,cutListEPICS); 
   if(rc!=0) return 1;

   // get vectors of data 
   std::vector<scalerData_t> raw,data; 
   mgr->GetVector_scaler("sbs",raw); 

   std::vector<epicsData_t> eraw,edata; 
   mgr->GetVector_epics(eraw); 
 
   std::string theVar;
   std::vector<double> out; 

   char var[200]; 
   producedVariable_t aPt;
   std::vector<producedVariable_t> pVar,pVarEPICS; 

   // apply cuts and compute stats  
   const int NC = cutList.size(); 
   for(int i=0;i<NC;i++){
      // define variable and get a vector of all data  
      sprintf(var,"%s.rate",devName.c_str()); 
      theVar = var; 
      // apply the cut 
      cut_util::ApplyCut(cutList[i],raw,data); 
      // compute stats on the cut data for the variable
      GetStats(theVar,data,out);
      // save results
      aPt.dev        = devName; 
      aPt.beam_state = cutList[i].beam_state; 
      aPt.group      = cutList[i].group; 
      aPt.time       = out[0]; 
      aPt.mean       = out[1]; 
      aPt.stdev      = out[2];
      pVar.push_back(aPt);  
      // set up for next cut 
      data.clear();
      out.clear();
   }

   // write results to file 
   char outpath[200];
   sprintf(outpath,"%s/%s.csv",out_dir,devName.c_str());
   bcm_util::WriteToFile(outpath,pVar);

   // apply cuts to EPICS data if necessary

   const int NCE = cutListEPICS.size(); 
   bool calcEPICS=false;

   // only do this for U1 and D1 signals
   std::string devNameEPICS="NONE";
   if(devName.compare("u1")==0){
      calcEPICS    = true; 
      devNameEPICS = "hac_bcm_dvm1_read"; 
   }else if(devName.compare("d1")==0){
      calcEPICS    = true; 
      devNameEPICS = "hac_bcm_dvm2_read"; 
   }

   if(calcEPICS){
      msg = "**** PROCESSING EPICS VARIABLE " + devNameEPICS + " ****";
      util_df::LogMessage(log_path,msg.c_str(),'a'); 
      for(int i=0;i<NCE;i++){
	 // define variable and get a vector of all data  
	 theVar = devNameEPICS; 
	 // apply the cut 
	 cut_util::ApplyCut(cutListEPICS[i],eraw,edata); 
	 // compute stats on the cut data for the variable
	 GetStats(theVar,edata,out);
	 // save results
	 aPt.dev        = devName; 
	 aPt.beam_state = cutListEPICS[i].beam_state; 
	 aPt.group      = cutListEPICS[i].group; 
	 aPt.time       = out[0]; 
	 aPt.mean       = out[1]; 
	 aPt.stdev      = out[2];
	 pVarEPICS.push_back(aPt);  
	 // set up for next cut 
	 edata.clear();
	 out.clear();
      }
      // write to file 
      sprintf(outpath,"%s/%s.csv",out_dir,devNameEPICS.c_str());
      bcm_util::WriteToFile(outpath,pVarEPICS);
   }
 
   return 0;
}
//______________________________________________________________________________
int GetStats(std::string var,std::vector<scalerData_t> data,std::vector<double> &out){
   // we want the time stamp, and mean BCM values for variable named var 
   
   const int N = data.size();
   std::vector<double> time,v; 
   for(int i=0;i<N;i++){
      time.push_back(data[i].time); 
      v.push_back( data[i].getValue(var) ); 
   }  

   // get values
   double theTime = time[0]; // note: using first time stamp! 
   double mean    = math_df::GetMean<double>(v);   
   double stdev   = math_df::GetStandardDeviation<double>(v);  

   // store in output vector
   out.push_back(theTime); 
   out.push_back(mean); 
   out.push_back(stdev); 
  
   return 0;
}
//______________________________________________________________________________
int GetStats(std::string var,std::vector<epicsData_t> data,std::vector<double> &out){
   // we want the time stamp, and mean BCM values for variable named var 
   
   const int N = data.size();
   std::vector<double> time,v; 
   for(int i=0;i<N;i++){
      time.push_back(data[i].time); 
      v.push_back( data[i].getValue(var) ); 
   }  

   // get values
   double theTime = time[0]; // note: using first time stamp! 
   double mean    = math_df::GetMean<double>(v);   
   double stdev   = math_df::GetStandardDeviation<double>(v);  

   // store in output vector
   out.push_back(theTime); 
   out.push_back(mean); 
   out.push_back(stdev); 
  
   return 0;
}

// Script for use with Unser calibration data
// Apply cuts to the data and produce an output file 

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"

#include "./include/codaRun.h"
#include "./src/BCMManager.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/Graph.cxx"

int unserProcessCuts(const char *runPath,const char *cutPath,const char *outpath){

   int rc=0;

   TString prefix;
   std::vector<codaRun_t> runList;  
   rc = bcm_util::LoadRuns(runPath,prefix,runList);
   if(rc!=0) return 1; 

   BCMManager *mgr = new BCMManager();

   TString filePath;  
   const int NR = runList.size();  
   for(int i=0;i<NR;i++){ 
      filePath = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d.root",
                      prefix.Data(),runList[i].runNumber,runList[i].stream,runList[i].segmentBegin,runList[i].segmentEnd);
      mgr->LoadFile(filePath,runList[i].runNumber);
   }

   std::vector<cut_t> cutList; 
   rc = bcm_util::LoadCuts(cutPath,cutList);
   if(rc!=0) return 1;

   double mean=0,stdev=0;
   std::vector<double> event,unser,time,ct;

   producedVariable_t pt; 
   std::vector<producedVariable_t> data;    

   int M=0;
   double theTime=0;
   // let's do cuts on each variable defined 
   int NC = cutList.size(); 
   for(int i=0;i<NC;i++){
      // define variable and get a vector of all data  
      mgr->GetVector(cutList[i].arm.c_str(),"event"     ,event);
      mgr->GetVector(cutList[i].arm.c_str(),"time"      ,time );
      mgr->GetVector(cutList[i].arm.c_str(),"unser.rate",unser);
      bcm_util::GetStatsWithCuts(event,unser,cutList[i].low,cutList[i].high,mean,stdev);
      bcm_util::ApplyCuts(cutList[i].low,cutList[i].high,event,time,ct); 
      std::cout << Form("[Cuts applied: cut lo = %.3lf, cut hi = %.3lf, group: %d]: %s mean = %.3lf, stdev = %.3lf",
	                cutList[i].low,cutList[i].high,cutList[i].group,"unser.rate",mean,stdev) << std::endl;
      // save results
      M = ct.size(); 
      // theTime       = math_df::GetMean<double>(ct); 
      theTime       = ct[M-1]; // use last time as the timestamp  
      pt.dev        = cutList[i].dev; 
      pt.beam_state = cutList[i].beam_state; 
      pt.group      = cutList[i].group;
      pt.time       = theTime;
      pt.mean       = mean; 
      pt.stdev      = stdev; 
      data.push_back(pt);  
      // set up for next cut 
      unser.clear();
      event.clear();
      ct.clear();
      time.clear();
   }

   bcm_util::WriteToFile(outpath,data);

   return 0;
}

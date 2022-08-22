// Plot all BCM data vs run number  

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"

#include "./include/bcm.h"
#include "./include/codaRun.h"
#include "./include/charge.h"
#include "./src/ABA.cxx"
#include "./src/CSVManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/JSONManager.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

int bcmCheckCharge(const char *confPath){

   // settings 
   bool logScale = false;

   gStyle->SetOptStat(0);

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix  = jmgr->GetValueFromKey_str("ROOTfile-path"); 
   std::string runPath = jmgr->GetValueFromKey_str("run-path");  
   
   BCMManager *mgr = new BCMManager("NONE","NONE",false);

   std::vector<codaRun_t> runList;  
   rc = util_df::LoadRunList(runPath.c_str(),prefix.c_str(),runList);
   if(rc!=0) return 1; 

   util_df::LoadBCMData(runList,mgr); 
   if(rc!=0) return 1; 

   std::vector<int> rr; 
   mgr->GetRunList(rr); 
   const int NNR = rr.size();

   const int N = 6; 
   TString var[N] = {"u1.current","unew.current","d1.current","d3.current","d10.current","dnew.current"};

   charge_t qData; 
   std::vector<scalerData_t> runData;

   for(int i=0;i<NNR;i++){
      // get data for the run
      mgr->GetVector_scaler("sbs",rr[i],runData);
      std::cout << "------------------------------------" << std::endl;
      std::cout << Form("Run %d: ",rr[i]) << std::endl;
      // loop over all BCMs and compute charge
      for(int j=0;j<N;j++){ 
	 bcm_util::GetCharge(var[j].Data(),runData,qData);
         std::cout << Form("   %s: totalTime = %.3lf sec (%.1lf min), Q = (%.3E ± %.3E) C",
                           var[j].Data(),qData.totalTime,qData.totalTime/60.,qData.value,qData.error) << std::endl;
      }
      // set up for next run
      runData.clear();
   }
   
   delete mgr;
   delete jmgr; 

   return 0;
}

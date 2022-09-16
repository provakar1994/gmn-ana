// Plot all BCM data vs run number  

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"

//#include "../../dflay/include/bcm.h"
#include "../../dflay/include/codaRun.h"
#include "../../dflay/include/charge.h"
#include "../../dflay/src/ABA.cxx"
#include "../../dflay/src/CSVManager.cxx"
#include "../../dflay/src/BCMManager.cxx"
#include "../../dflay/src/Utilities.cxx"
#include "../../dflay/src/JSONManager.cxx"
#include "../../dflay/src/Graph.cxx"
#include "../../dflay/src/bcmUtilities.cxx"
#include "../../dflay/src/cutUtilities.cxx"

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

   charge_t qData; 
   std::vector<scalerData_t> runData;

   TCanvas *c1 = new TCanvas("c1", "c1", 600, 800);
   c1->Divide(1,2);
   c1->cd(1);

   for(int i=0;i<NNR;i++){

     // get data for the run
     mgr->GetVector_scaler("sbs",rr[i],runData);

     if (i == 0) {
       TGraph *g1 = bcm_util::GetTGraph("event", "dnew.current", runData);
       graph_df::SetLabels(g1, Form("Run # %d: Beam Current (uA)", runData[0].runNumber), 
			   "event", "dnew.current (uA)");
       g1->Draw("alp");

       // c1->cd(2);
       // TGraph *g2 = bcm_util::GetTGraph_timeStep("event", runData);
       // graph_df::SetLabels(g2, Form("Run # %d: Time steps", runData[0].runNumber), 
       // 			   "event", "Time steps (s)");
       // g2->Draw("alp");
       // TH1D *h1 = bcm_util::GetTH1D(runData, "BBCalHi.scalerRate", 500, 3500, 5500);
       // h1->Draw();

       // c1->cd(2); 
       // // TGraph *g4 = bcm_util::GetTGraph_charge("event", "dnew.current", runData); // dflay method
       // TGraph *g4 = bcm_util::GetTGraph_charge_pd("dnew.cnt", "event", runData);
       // graph_df::SetLabels(g4, Form("Run # %d: Beam Charge (C) using dnew source", runData[0].runNumber),
       // 			   "event", "Charge (C)");
       // g4->Draw("alp");

       // c1->cd(2);
       // TGraph *g3 = bcm_util::GetTGraph_charge_pd("dnew.cnt", "time103kHz", runData);
       // graph_df::SetLabels(g3, Form("Run # %d: Beam Charge (C) using dnew source", runData[0].runNumber), 
       // 			   "time (s)", "Charge (C)");
       // g3->Draw("alp");

       // c1->cd(2);
       // TGraph *g4 = bcm_util::GetTGraph_charge_pd("dnew.cnt", "event", runData);
       // graph_df::SetLabels(g4, "Beam Charge (C) using dnew source", "event", "Charge (C)");
       // g4->Draw("alp");

       c1->cd(2);
       TH1D *h1 = bcm_util::GetTH1D(runData, "liveTime", 100, 0.5, 1.5);
       // std::cout << " Live Time " << 0.5 + h1->GetMaximumBin()*h1->GetBinWidth(1) << std::endl;
       h1->Draw();

     }

     // c1->cd(3);
     // TGraph *g4 = bcm_util::GetTGraph("event", "L1A.scalerRate", runData);
     // //graph_df::SetLabels(g3, "Beam Charge (C) using dnew.current", "event", "Charge (C)");
     // g4->Draw("alp");

     // c1->cd(3);
     // TH1D *h1 = bcm_util::GetTH1D(runData, "BBCalHi.scalerRate", 500, 3500, 5500);
     // //graph_df::SetLabels(g3, "Beam Charge (C) using dnew.current", "event", "Charge (C)");
     // //h1->Draw();
 
     // TH2D *h2 = bcm_util::GetTH2D(runData, "time", "dnew.current", 1000, 0, 1000, 10, 0, 5);
     // h2->Draw("colz");

     std::cout << "------------------------------------" << std::endl;
     std::cout << Form("Run %d: ",rr[i]) << std::endl;

     bcm_util::GetCharge_pd("dnew.cnt", runData, qData);
     std::cout << Form(" dnew.cnt : totalTime = %.3lf sec (%.1lf min), Q = %.3E C",
		       qData.totalTime,qData.totalTime/60.,qData.value) << std::endl;
     std::cout << " DAQ Live Time = " << bcm_util::GetDAQLiveTime(runData) << std::endl;

     // set up for next run
     runData.clear();
   }

   delete mgr;
   delete jmgr; 

   return 0;
}

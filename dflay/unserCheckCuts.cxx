// Script for use with Unser calibration data
// Plot the Unser data and the cuts 

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

int unserCheckCuts(const char *runPath){

   // settings 
   bool logScale   = false;

   gStyle->SetOptStat(0);

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

   // mgr->Print("sbs"); 

   std::vector<cut_t> cutList; 
   rc = bcm_util::LoadCuts("./input/cut-list_unser-calib.csv",cutList);
   if(rc!=0) return 1;

   double mean=0,stdev=0;
   std::vector<double> event,unser;

   const int NC = cutList.size();
   TLine **lo = new TLine*[NC];
   TLine **hi = new TLine*[NC];

   double yMin = 770E+3; 
   double yMax = 1200E+3; 

   // let's do cuts on each variable defined 
   for(int i=0;i<NC;i++){
      // define variable and get a vector of all data  
      mgr->GetVector(cutList[i].arm.c_str(),"event"     ,event);
      mgr->GetVector(cutList[i].arm.c_str(),"unser.rate",unser);
      bcm_util::GetStatsWithCuts(event,unser,cutList[i].low,cutList[i].high,mean,stdev);
      std::cout << Form("[Cuts applied: cut lo = %.3lf, cut hi = %.3lf, group: %d]: %s mean = %.3lf, stdev = %.3lf",
	    cutList[i].low,cutList[i].high,cutList[i].group,"unser.rate",mean,stdev) << std::endl;
      // make lines we can plot 
      lo[i] = new TLine(cutList[i].low ,yMin,cutList[i].low,yMax);
      lo[i]->SetLineWidth(2); 
      hi[i] = new TLine(cutList[i].high,yMin,cutList[i].high,yMax);
      hi[i]->SetLineWidth(2); 
      if(cutList[i].beam_state.compare("on")==0){
	 lo[i]->SetLineColor(kGreen+1); 
	 hi[i]->SetLineColor(kGreen+1); 
      }else{
	 lo[i]->SetLineColor(kRed); 
	 hi[i]->SetLineColor(kRed); 
      }

      // set up for next cut 
      unser.clear();
      event.clear();
   }

   TGraph *ge = mgr->GetTGraph("sbs","event","unser.rate");
   graph_df::SetParameters(ge,20,kBlack); 
 
   TGraph *gt = mgr->GetTGraph("sbs","time","unser.rate");
   graph_df::SetParameters(gt,20,kBlack); 

   TString Title      = Form("Unser Calibration Data");
   TString yAxisTitle = Form("Unser Rate [Hz]");

   TCanvas *c1 = new TCanvas("c1","Unser Calibration",1200,800);

   c1->cd();
   ge->Draw("alp");
   graph_df::SetLabels(ge,Title,"event",yAxisTitle);
   ge->Draw("alp");
   for(int i=0;i<NC;i++){
      lo[i]->Draw("same"); 
      hi[i]->Draw("same"); 
   } 
   c1->Update();

   // TCanvas *c2 = new TCanvas("c2","Unser Calibration Data vs Time",1200,800); 

   // c2->cd();
   // gt->Draw("alp");
   // graph_df::SetLabels(gt,Title,"Time [sec]",yAxisTitle); 
   // gt->Draw("alp");
   // c2->Update();

   return 0;
}

// Script for use with Unser calibration data
// Plot the Unser data  

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

int unserPlot(const char *runPath){

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
   c1->Update();

   TCanvas *c2 = new TCanvas("c2","Unser Calibration Data vs Time",1200,800); 

   c2->cd();
   gt->Draw("alp");
   graph_df::SetLabels(gt,Title,"Time [sec]",yAxisTitle); 
   gt->Draw("alp");
   c2->Update();

   return 0;
}

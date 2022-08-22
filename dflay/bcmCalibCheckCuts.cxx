// Check the chosen cuts for the BCM data  

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
#include "./src/CSVManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/JSONManager.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

int bcmCalibCheckCuts(const char *confPath){

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix  = jmgr->GetValueFromKey_str("ROOTfile-path");
   std::string tag     = jmgr->GetValueFromKey_str("tag");
   delete jmgr;

   char run_path[200],cut_path[200],out_dir[200];
   sprintf(run_path,"./input/%s/runlist.csv"     ,tag.c_str());
   sprintf(cut_path,"./input/%s/cuts/cutlist.csv",tag.c_str());
   sprintf(out_dir ,"./output/%s"                ,tag.c_str());

   BCMManager *mgr = new BCMManager("NONE","NONE",false);

   std::vector<codaRun_t> runList;
   rc = util_df::LoadRunList(run_path,prefix.c_str(),runList);
   if(rc!=0) return 1;

   util_df::LoadBCMData(runList,mgr);
   if(rc!=0) return 1;

   std::vector<cut_t> cutList;
   rc = bcm_util::LoadCuts(cut_path,cutList);
   if(rc!=0) return 1;

   const int N = 7; 
   TString varName[N] = {"unser.rate","u1.rate","unew.rate","dnew.rate","d1.rate","d3.rate","d10.rate"};

   TString xAxis = "runEvent"; 

   TString theVar,theCutVar; 
   double mean=0,stdev=0;
   std::vector<double> time,v; 

   const int NC = cutList.size();
   TLine **lo = new TLine*[NC]; 
   TLine **hi = new TLine*[NC]; 

   int color=0;
  
   // let's do cuts on each variable defined 
   for(int i=0;i<NC;i++){
      // define variable and get a vector of all data 
      theCutVar = Form("%s",cutList[i].cut_var.c_str());
      if(cutList[i].arm.compare("E")==0){
	 // EPICS variable 
	 theVar = Form("%s",cutList[i].dev.c_str());
      }else{
	 theVar = Form("%s.rate",cutList[i].dev.c_str());
      }
      mgr->GetVector(cutList[i].arm.c_str(),theCutVar.Data(),time); 
      mgr->GetVector(cutList[i].arm.c_str(),theVar.Data()   ,v); 
      bcm_util::GetStatsWithCuts(time,v,cutList[i].low,cutList[i].high,mean,stdev);
      std::cout << Form("[Cuts applied: cut lo = %.3lf, cut hi = %.3lf, group: %d]: %s mean = %.3lf, stdev = %.3lf",
                        cutList[i].low,cutList[i].high,cutList[i].group,theVar.Data(),mean,stdev) << std::endl; 
      // make lines we can plot
      color = kRed;
      if(i%2==0) color = kGreen+2; 
      lo[i] = new TLine(cutList[i].low ,0,cutList[i].low ,900E+3);
      lo[i]->SetLineColor(color); 
      hi[i] = new TLine(cutList[i].high,0,cutList[i].high,900E+3);
      hi[i]->SetLineColor(color);
      // set up for next cut 
      v.clear();
      time.clear();
   }
 
   // create histos and TGraphs 
   TGraph **g = new TGraph*[N];
   TCanvas **c = new TCanvas*[N]; 

   for(int i=0;i<N;i++){
      g[i] = mgr->GetTGraph("sbs",xAxis,varName[i]);
      graph_df::SetParameters(g[i],20,kBlack); 
   }

   TString cName,cTitle;

   for(int i=0;i<N;i++){
      cName  = Form("c%d",i);
      cTitle = Form("%s",varName[i].Data());
      c[i]   = new TCanvas(cName,cTitle,1200,800);
      c[i]->cd();
      g[i]->Draw("ap");
      graph_df::SetLabels(g[i],varName[i].Data(),xAxis.Data(),varName[i].Data());
      g[i]->Draw("ap");
      for(int j=0;j<NC;j++){
	 lo[j]->Draw("same"); 
	 hi[j]->Draw("same"); 
      }
      c[i]->Update();
   }

   TGraph *gEPICSCurrent = mgr->GetTGraph("E","event","IBC1H04CRCUR2"); 
   gEPICSCurrent->SetMarkerStyle(20);

   // EPICS plots

   TCanvas *ce = new TCanvas("ce","EPICS Beam Current",1200,800); 

   ce->cd();
   gEPICSCurrent->Draw("ap");
   graph_df::SetLabels(gEPICSCurrent,"","event","IBC1H04CRCUR2 [#muA]"); 
   gEPICSCurrent->Draw("ap");
   ce->Update();  

   return 0;
}

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
#include "./src/ABA.cxx"
#include "./src/CSVManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/JSONManager.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

#define MICROAMPS 1E-6

TString GetTitle(std::vector<int> run); 
int checkEvt(std::vector<scalerData_t> data); 

int bcmCheckCurrent(const char *confPath){

   // settings 
   bool logScale = false;

   gStyle->SetOptStat(0);

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix       = jmgr->GetValueFromKey_str("ROOTfile-path"); 
   std::string runPath      = jmgr->GetValueFromKey_str("run-path");  
   std::string cutPath      = jmgr->GetValueFromKey_str("cut-path"); 
   
   BCMManager *mgr = new BCMManager("NONE","NONE",false);

   std::vector<codaRun_t> runList;  
   rc = util_df::LoadRunList(runPath.c_str(),prefix.c_str(),runList);
   if(rc!=0) return 1; 

   util_df::LoadBCMData(runList,mgr); 
   if(rc!=0) return 1; 

   // do stats by run number 
   std::vector<scalerData_t> sbs,lhrs;
   mgr->GetVector_scaler("sbs" ,sbs);  
   mgr->GetVector_scaler("Left",lhrs);  

   std::vector<epicsData_t> edata; 
   mgr->GetVector_epics(edata);  

   // // load cuts 
   // std::vector<cut_t> cutList; 
   // cut_util::LoadCuts(cutPath.c_str(),cutList);
   // // apply cuts  
   // cut_util::ApplyCuts(cutList,rawData,data);

   // sort the data by run number 
   std::sort(sbs.begin() ,sbs.end() ,compareScalerData_byRun); 
   std::sort(lhrs.begin(),lhrs.end(),compareScalerData_byRun); 
   std::sort(edata.begin(),edata.end(),compareEPICSData_byRun); 

   checkEvt(sbs); 

   std::vector<int> rr; 
   mgr->GetRunList(rr); 
   const int NNR = rr.size();

   const int N = 6; 
   TString var[N]  = {"u1.cnt","unew.cnt","d1.cnt","d3.cnt","d10.cnt","dnew.cnt"};
   TString varR[N] = {"u1.rate","unew.rate","d1.rate","d3.rate","d10.rate","dnew.rate"};
   TString varC[N] = {"u1.current","unew.current","d1.current","d3.current","d10.current","dnew.current"};
   int color[N]    = {kBlack,kRed,kBlue,kGreen+2,kMagenta,kCyan+2}; 

   std::string xVar    = "runEvent"; 
   TString xAxisTitle  = Form("");
   TString canvasTitle = GetTitle(rr);
   TString Title,yAxisTitle; 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   TMultiGraph **mgCur = new TMultiGraph*[N]; 
   TMultiGraph **mgCnt = new TMultiGraph*[N]; 
   TMultiGraph **mgRte = new TMultiGraph*[N]; 
   TGraph **gLeft_cur  = new TGraph*[N]; 
   TGraph **gSBS_cur   = new TGraph*[N]; 
   TGraph **gLeft_cnt  = new TGraph*[N]; 
   TGraph **gSBS_cnt   = new TGraph*[N]; 
   TGraph **gLeft_rte  = new TGraph*[N]; 
   TGraph **gSBS_rte   = new TGraph*[N]; 

   TCanvas *c1 = new TCanvas("c1","BCM Currents",1200,800);
   c1->Divide(3,2);
 
   TCanvas *c2 = new TCanvas("c2","BCM Counts",1200,800);
   c2->Divide(3,2);  

   TCanvas *c3 = new TCanvas("c3","BCM Rates" ,1200,800);
   c3->Divide(3,2); 

   // get BCM plots 
   for(int i=0;i<N;i++){
      // currents
      Title      = Form("%s",varC[i].Data());
      yAxisTitle = Form("%s [#muA]",varC[i].Data());
      gSBS_cur[i]  = bcm_util::GetTGraph(xVar.c_str(),varC[i].Data(),sbs ); 
      gLeft_cur[i] = bcm_util::GetTGraph(xVar.c_str(),varC[i].Data(),lhrs);
      graph_df::SetParameters(gSBS_cur[i] ,21,kBlack);
      graph_df::SetParameters(gLeft_cur[i],20,kRed  );
      mgCur[i] = new TMultiGraph();
      mgCur[i]->Add(gSBS_cur[i],"p");  
      mgCur[i]->Add(gLeft_cur[i],"p");
      if(i==0){
	 L->AddEntry(gSBS_cur[i] ,"SBS" ,"p"); 
	 L->AddEntry(gLeft_cur[i],"LHRS","p"); 
      } 
      c1->cd(i+1);
      mgCur[i]->Draw("a"); 
      graph_df::SetLabels(mgCur[i],Title,xAxisTitle,yAxisTitle); 
      if(xVar.compare("time")==0) graph_df::UseTimeDisplay(mgCur[i]); 
      mgCur[i]->Draw("a");
      if(i==0) L->Draw("same"); 
      c1->Update();
      // counts
      Title      = Form("%s",var[i].Data());
      yAxisTitle = Form("%s",var[i].Data());
      gSBS_cnt[i]  = bcm_util::GetTGraph(xVar.c_str(),var[i].Data(),sbs ); 
      gLeft_cnt[i] = bcm_util::GetTGraph(xVar.c_str(),var[i].Data(),lhrs);
      graph_df::SetParameters(gSBS_cnt[i] ,21,kBlack);
      graph_df::SetParameters(gLeft_cnt[i],20,kRed  );
      mgCnt[i] = new TMultiGraph();
      mgCnt[i]->Add(gSBS_cnt[i],"p");  
      mgCnt[i]->Add(gLeft_cnt[i],"p");
      c2->cd(i+1);
      mgCnt[i]->Draw("a"); 
      graph_df::SetLabels(mgCnt[i],Title,xAxisTitle,yAxisTitle); 
      if(xVar.compare("time")==0) graph_df::UseTimeDisplay(mgCnt[i]); 
      mgCnt[i]->Draw("a");
      if(i==0) L->Draw("same"); 
      c2->Update();
      // rates
      Title      = Form("%s",varR[i].Data());
      yAxisTitle = Form("%s [Hz]",varR[i].Data());
      gSBS_rte[i]  = bcm_util::GetTGraph(xVar.c_str(),varR[i].Data(),sbs ); 
      gLeft_rte[i] = bcm_util::GetTGraph(xVar.c_str(),varR[i].Data(),lhrs);
      graph_df::SetParameters(gSBS_rte[i] ,21,kBlack);
      graph_df::SetParameters(gLeft_rte[i],20,kRed  );
      mgRte[i] = new TMultiGraph();
      mgRte[i]->Add(gSBS_rte[i],"p");  
      mgRte[i]->Add(gLeft_rte[i],"p");
      c3->cd(i+1);
      mgRte[i]->Draw("a"); 
      graph_df::SetLabels(mgRte[i],Title,xAxisTitle,yAxisTitle); 
      if(xVar.compare("time")==0) graph_df::UseTimeDisplay(mgRte[i]); 
      mgRte[i]->Draw("a");
      if(i==0) L->Draw("same"); 
      c3->Update();
   }
  
   // TCanvas *c3 = new TCanvas("c3",canvasTitle,1200,800);

   // c3->cd();
   // mgc->Draw("a"); 
   // graph_df::SetLabels(mgc,Title,xAxisTitle,"Beam Current [#muA]");
   // // graph_df::SetLabelSizes(mgc,0.05,0.06);  
   // graph_df::UseTimeDisplay(mgc);  
   // mgc->Draw("a");
   // L->Draw("same"); 
   // c3->Update(); 

   delete mgr;
   delete jmgr; 

   return 0;
}
//______________________________________________________________________________
int checkEvt(std::vector<scalerData_t> data){
   std::cout << "Checking trigger event..." << std::endl;
   std::vector<double> ev1,ev2;
   const int N = data.size();
   for(int i=0;i<N;i++){
      if(data[i].triggerEvent!=data[i].triggerEvent2){
	 std::cout << Form("DISCREPANCY: scaler evt = %d, triggerEvt = %lld, triggerEvt2 = %lld",
                           i,data[i].triggerEvent,data[i].triggerEvent2) << std::endl;
      }
      ev1.push_back( (double)data[i].triggerEvent  ); 
      ev2.push_back( (double)data[i].triggerEvent2 ); 
   }
   std::cout << "--> Done." << std::endl;

   // make a plot
   TGraph *g = graph_df::GetTGraph(ev1,ev2); 
   graph_df::SetParameters(g,20,kBlack);

   TCanvas *cc = new TCanvas("cc","Event Check",1000,800);

   cc->cd();
   g->Draw("alp"); 
   graph_df::SetLabels(g,"Trigger Event Check","evnum","evNumber");  
   g->Draw("alp"); 
   cc->Update();

   return 0; 
}
//______________________________________________________________________________
TString GetTitle(std::vector<int> run){
   TString title="";
   int first=0,last=0;
   const int N = run.size();
   if(N==0){
      std::cout << "[GetTitle]: No runs!" << std::endl;
   }else if(N==1){
      title = Form("Run %d: Beam Current",run[0]);
   }else{
      title = Form("Run %d--%d: Beam Current",run[0],run[N-1]);
   } 
   return title;
}

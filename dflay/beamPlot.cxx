// Plot all BCM data  

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"

#include "./include/scalerData.h"
#include "./include/epicsData.h"
#include "./include/producedVariable.h"
#include "./include/calibCoeff.h"
#include "./include/codaRun.h"
#include "./src/ABA.cxx"
#include "./src/Graph.cxx"
#include "./src/CSVManager.cxx"
#include "./src/JSONManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/BeamManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/bcmUtilities.cxx"

int beamPlot(const char *confPath){

   // settings 
   bool logScale   = false;

   gStyle->SetOptStat(0);

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix  = jmgr->GetValueFromKey_str("ROOTfile-path");
   std::string runPath = jmgr->GetValueFromKey_str("run-path");
   delete jmgr; 

   std::vector<codaRun_t> runList;  
   rc = bcm_util::LoadRuns(runPath.c_str(),runList);
   if(rc!=0) return 1; 

   BeamManager *mgr = new BeamManager();

   std::vector<int> md; 

   TString filePath;  
   const int NR = runList.size();  
   for(int i=0;i<NR;i++){ 
      util_df::GetROOTFileMetaData(prefix.c_str(),runList[i].runNumber,md);
      runList[i].stream       = md[0];
      runList[i].segmentBegin = md[1];
      runList[i].segmentEnd   = md[2];
      filePath = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d.root",
                      prefix.c_str(),runList[i].runNumber,runList[i].stream,runList[i].segmentBegin,runList[i].segmentEnd);
      mgr->LoadFile(filePath,runList[i].runNumber);
      md.clear();
   }

   const int N = 8; 
   TString varName[N] = {"Raster.rawcur.x","Raster.rawcur.y","Raster2.rawcur.x","Raster2.rawcur.y",
                         "BPMA.rawcur.1"  ,"BPMA.rawcur.2"  ,"BPMB.rawcur.1"   ,"BPMB.rawcur.2"};
   TString xAxis      = Form("event"); 

   // create graphs and histograms 

   TGraph **g = new TGraph*[N];
   TH1F **h   = new TH1F*[N]; 

   int NBin = 100; 
   // histo bounds  u1  unew    unser  dnew d1     d3    d10  
   double min[N] = {0     ,0     ,0     ,0     ,0     ,0     ,0     };
   double max[N] = {100E+3,100E+3,100E+3,100E+3,100E+3,100E+3,100E+3};

   for(int i=0;i<N;i++){
      g[i] = mgr->GetTGraph("SBSrb",xAxis.Data(),varName[i]);
      graph_df::SetParameters(g[i],20,kBlack);
      h[i] = mgr->GetTH1F("SBSrb",varName[i].Data(),NBin,min[i],max[i]); 
   }

   // 2D raster plot, BPM plot
   int NBin2D = 4000;  
   TH2F *r2D_1 = mgr->GetTH2F("SBSrb","Raster1.rawcur.x"   ,"Raster1.rawcur.y" ,NBin2D,0,95E+3,NBin2D,0,95E+3);
   TH2F *r2D_2 = mgr->GetTH2F("SBSrb","Raster2.rawcur.x"   ,"Raster2.rawcur.y",NBin2D,0,95E+3,NBin2D,0,95E+3);
   TH2F *p2D_A = mgr->GetTH2F("SBSrb","BPMA.x"  ,"BPMA.y"  ,NBin2D,-10E-3,10E-3,NBin2D,-10E-3,10E-3);
   TH2F *p2D_B = mgr->GetTH2F("SBSrb","BPMB.x"  ,"BPMB.y"  ,NBin2D,-10E-3,10E-3,NBin2D,-10E-3,10E-3);
   TH2F *tgt   = mgr->GetTH2F("SBSrb","target.x","target.y",NBin2D,-4,4,NBin2D,-4,4);

   // epics data
   const int NE = 4; 
   TGraph **ge = new TGraph*[NE];  
   TH1F **he   = new TH1F*[NE];  
   TString epicsVarName[NE] = {"IPM1H04A_XPOS","IPM1H04A_YPOS","IPM1H04E_XPOS","IPM1H04E_YPOS"};
   for(int i=0;i<NE;i++){
      ge[i] = mgr->GetTGraph("E",xAxis.Data(),epicsVarName[i].Data()); 
      graph_df::SetParameters(ge[i],20,kBlack);
      he[i] = mgr->GetTH1F("E",epicsVarName[i].Data(),NBin,min[i],max[i]); 
   } 

   TCanvas *c1a = new TCanvas("c1a","Beam Check",1200,800);
   c1a->Divide(2,2);

   TCanvas *c1b = new TCanvas("c1b","Beam Check",1200,800);
   c1b->Divide(2,2);

   TString Title,xAxisTitle,yAxisTitle;
   if(xAxis=="time"){
      xAxisTitle = Form("%s [sec]",xAxis.Data());
   }else{
      xAxisTitle = Form("%s",xAxis.Data());
   }

   for(int i=0;i<N/2;i++){
      c1a->cd(i+1);
      Title      = Form("%s"     ,varName[i].Data());
      yAxisTitle = Form("%s [Hz]",varName[i].Data());
      g[i]->Draw("alp");
      graph_df::SetLabels(g[i],Title,xAxisTitle,yAxisTitle); 
      g[i]->Draw("alp");
      c1a->Update();
      // next canvas 
      c1b->cd(i+1);
      Title      = Form("%s"     ,varName[i+3].Data());
      yAxisTitle = Form("%s [Hz]",varName[i+3].Data());
      g[i+3]->Draw("alp");
      graph_df::SetLabels(g[i+3],Title,xAxisTitle,yAxisTitle); 
      g[i+3]->Draw("alp");
      c1b->Update();
   }

   // last one 
   Title      = Form("%s"     ,varName[6].Data());
   yAxisTitle = Form("%s [Hz]",varName[6].Data());
   c1b->cd(4); 
   g[6]->Draw("alp");
   graph_df::SetLabels(g[6],Title,xAxisTitle,yAxisTitle); 
   g[6]->Draw("alp");
   c1b->Update();

   // histograms
   TCanvas *c2a = new TCanvas("c2a","Beam Check [Histograms]",1200,800);
   c2a->Divide(2,2);

   TCanvas *c2b = new TCanvas("c2b","Beam Check [Histograms]",1200,800);
   c2b->Divide(2,2);

   for(int i=0;i<N/2;i++){
      c2a->cd(i+1);
      h[i]->Draw("");
      c2a->Update();
      // next canvas 
      c2b->cd(i+1);
      h[i+3]->Draw("");
      c2b->Update();
   }

   // last one 
   c2b->cd(4); 
   h[6]->Draw("");
   c2b->Update();

   TCanvas *c3 = new TCanvas("c3","Raster and BPM Plots",1200,800);
   c3->Divide(2,2);

   c3->cd(1); 
   r2D_1->Draw("colz");
   c3->Update();

   c3->cd(2); 
   r2D_2->Draw("colz");
   c3->Update();

   c3->cd(3); 
   p2D_A->Draw("colz");
   c3->Update();

   c3->cd(4); 
   p2D_B->Draw("colz");
   c3->Update();

   TCanvas *c4 = new TCanvas("c4","Beam Position at Target",1200,800);

   c4->cd(); 
   tgt->Draw("colz");
   c4->Update();

   TCanvas *c5 = new TCanvas("c5","EPICS BPM Data",1200,800);
   c5->Divide(2,2); 

   for(int i=0;i<NE;i++){
      c5->cd(i+1); 
      ge[i]->Draw("ap");
      graph_df::SetLabels(ge[i],epicsVarName[i],xAxis,epicsVarName[i]); 
      ge[i]->Draw("ap");
      c5->Update();
   }

   TCanvas *c6 = new TCanvas("c6","EPICS BPM Data Histos",1200,800);
   c6->Divide(2,2); 

   for(int i=0;i<NE;i++){
      c6->cd(i+1); 
      he[i]->Draw("");
      // graph_df::SetLabels(ge[i],epicsVarName[i],xAxis,epicsVarName[i]); 
      he[i]->Draw("");
      c6->Update();
   }



   return 0;
}

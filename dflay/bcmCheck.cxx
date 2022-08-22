// testing the BCM data 

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"

bool gIsDebug = false;

int GetEvents(int run,std::string system,TTree *T); 

int bcmCheck(){

   // settings 
   bool logScale   = false;
   bool sameCanvas = true;
   
   gStyle->SetOptStat(0);
   
   int run;
   std::cout << "Enter run number: ";
   std::cin  >> run;

   int startSegment = 0; 
   int endSegment   = 0; 

   TString prefix   = Form("/lustre19/expphy/volatile/halla/sbs/flay/GMnAnalysis/rootfiles");
   TString filePath = Form("%s/gmn_replayed-beam_%d_stream0_seg%d_%d.root",prefix.Data(),run,startSegment,endSegment);

   TString treePath_T    = Form("%s/T"     ,filePath.Data()); 
   TString treePath_LHRS = Form("%s/TSLeft",filePath.Data()); 
   TString treePath_SBS  = Form("%s/TSsbs" ,filePath.Data()); 

   int NEntries=0,NTrees=0;

   TChain *chT = new TChain("T"); 
   chT->Add(treePath_T);
   NEntries = chT->GetEntries();
   TTree *T = chT->GetTree();
 
   TChain *chL = new TChain("TSLeft"); 
   chL->Add(treePath_LHRS); 
   NEntries = chL->GetEntries();
   TTree *TSLeft = chL->GetTree();

   TChain *chS = new TChain("TSsbs"); 
   chS->Add(treePath_SBS); 
   NEntries = chS->GetEntries();
   TTree *TSsbs = chS->GetTree();
   
   if(gIsDebug) std::cout << "Tree addresses: T = " << T << " TSLeft = " << TSLeft << " TSsbs = " << TSsbs << std::endl;

   const int N = 8;

   int NBin = 100; 

   // histo bounds        u1  unew unser dnew d1 d3 d10 103.7 kHz clk
   double loBeamOff[N] = {50  ,0  ,700E+3,0    ,1E+3  ,0  ,0  ,102E+3}; 
   double hiBeamOff[N] = {150,100 ,850E+3,20E+3,1.5E+3,100,100,105E+3}; 

   TString leftVarRate[N] = {"Left.bcm.u1.rate","Left.bcm.unew.rate","Left.bcm.unser.rate",
                             "Left.bcm.dnew.rate","Left.bcm.d1.rate","Left.bcm.d3.rate","Left.bcm.d10.rate",
                             "Left.104kHz_CLK.rate"};  

   TString leftVarCnt[N]  = {"Left.bcm.u1.cnt","Left.bcm.unew.cnt","Left.bcm.unser.cnt",
                             "Left.bcm.dnew.cnt","Left.bcm.d1.cnt","Left.bcm.d3.cnt","Left.bcm.d10.cnt",  
                             "Left.104kHz_CLK.cnt"};  

   TString sbsVarRate[N]  = {"sbs.bcm.u1.rate"  ,"sbs.bcm.unew.rate","sbs.bcm.unser.rate",
                             "sbs.bcm.dnew.rate","sbs.bcm.d1.rate"  ,"sbs.bcm.d3.rate","sbs.bcm.d10.rate",
                             "sbs.104kHz_CLK.rate"};  

   TString sbsVarCnt[N]   = {"sbs.bcm.u1.cnt"  ,"sbs.bcm.unew.cnt","sbs.bcm.unser.cnt",
                             "sbs.bcm.dnew.cnt","sbs.bcm.d1.cnt"  ,"sbs.bcm.d3.cnt","sbs.bcm.d10.cnt", 
                             "sbs.104kHz_CLK.cnt"};  

   TString Title[N]       = {"u1 Rate [Hz]"  ,"unew Rate [Hz]","Unser Rate [Hz]",
                             "dnew Rate [Hz]","d1 Rate [Hz]"  ,"d3 Rate [Hz]"   ,"d10 Rate [Hz]","103.7 kHz Clock [Hz]"};

   double timeMax = 150;  

   TString leftClkVar = Form("Left.104kHz_CLK.cnt/Left.104kHz_CLK.rate");
   TString sbsClkVar  = Form("sbs.BBCalHi.RF.scaler/sbs.BBCalHi.RF.scalerRate");

   TString histoStr,varStr; 
   for(int i=0;i<N;i++){
      // LHRS
      // rates
      histoStr = Form("h%d(%d,%.0lf,%.0lf)",i+1,NBin,loBeamOff[i],hiBeamOff[i]);
      varStr   = leftVarRate[i]; 
      TSLeft->Project(histoStr,varStr,"","");
      // counts 
      histoStr = Form("h%dc",i+1); 
      varStr   = leftVarCnt[i]; 
      TSLeft->Project(histoStr,varStr,"","");
      // rate vs time  
      histoStr = Form("h%dvt(%d,0,%.0lf,%d,%.0lf,%.0lf)",i+1,NBin,timeMax,NBin,loBeamOff[i],hiBeamOff[i]); 
      varStr   = Form("%s:%s",leftVarRate[i].Data(),leftClkVar.Data()); 
      TSLeft->Project(histoStr,varStr,"","");
      // SBS
      // rates
      histoStr = Form("g%d(%d,%.0lf,%.0lf)",i+1,NBin,loBeamOff[i],hiBeamOff[i]);
      varStr   = sbsVarRate[i]; 
      TSsbs->Project(histoStr,varStr,"","");
      // counts  
      histoStr = Form("g%dc",i+1); 
      varStr   = sbsVarCnt[i]; 
      TSsbs->Project(histoStr,varStr,"","");
      // rate vs time  
      histoStr = Form("g%dvt(%d,0,%.01lf,%d,%.0lf,%.0lf)",i+1,NBin,timeMax,NBin,loBeamOff[i],hiBeamOff[i]); 
      varStr   = Form("%s:%s",sbsVarRate[i].Data(),sbsClkVar.Data()); 
      TSsbs->Project(histoStr,varStr,"","");
   }

   TH1F **h   = new TH1F*[N];
   TH1F **hc  = new TH1F*[N];
   TH2F **hvt = new TH2F*[N];
   TH1F **g   = new TH1F*[N];
   TH1F **gc  = new TH1F*[N];
   TH2F **gvt = new TH2F*[N];

   for(int i=0;i<N;i++){
      // LHRS 
      h[i] = (TH1F *)gDirectory->Get( Form("h%d",i+1) );
      h[i]->SetMarkerColor(kBlue);
      h[i]->SetMarkerStyle(20);
      h[i]->SetLineColor(kBlue);
      h[i]->SetLineWidth(1);
      h[i]->SetTitle(Title[i]); 
      hc[i] = (TH1F *)gDirectory->Get( Form("h%dc",i+1) );
      hc[i]->SetMarkerColor(kBlue);
      hc[i]->SetMarkerStyle(20);
      hc[i]->SetLineColor(kBlue);
      hc[i]->SetLineWidth(1);
      hvt[i] = (TH2F *)gDirectory->Get( Form("h%dvt",i+1) );
      hvt[i]->SetMarkerColor(kBlue);
      hvt[i]->SetMarkerStyle(20);
      hvt[i]->SetLineColor(kBlue);
      hvt[i]->SetLineWidth(1);
      hvt[i]->SetTitle(Title[i]); 
      hvt[i]->GetXaxis()->SetTitle("Time [sec]"); 
      hvt[i]->GetXaxis()->CenterTitle(); 
      hvt[i]->GetYaxis()->SetTitle(Title[i]); 
      hvt[i]->GetYaxis()->CenterTitle(); 
      // SBS 
      g[i] = (TH1F *)gDirectory->Get( Form("g%d",i+1) );
      g[i]->SetMarkerColor(kRed);
      g[i]->SetMarkerStyle(20);
      g[i]->SetLineColor(kRed);
      g[i]->SetLineWidth(1);
      g[i]->SetTitle(Title[i]); 
      gc[i] = (TH1F *)gDirectory->Get( Form("g%dc",i+1) );
      gc[i]->SetMarkerColor(kRed);
      gc[i]->SetMarkerStyle(20);
      gc[i]->SetLineColor(kRed);
      gc[i]->SetLineWidth(1);
      gvt[i] = (TH2F *)gDirectory->Get( Form("g%dvt",i+1) );
      gvt[i]->SetMarkerColor(kRed);
      gvt[i]->SetMarkerStyle(20);
      gvt[i]->SetLineColor(kRed);
      gvt[i]->SetLineWidth(1);
      gvt[i]->SetTitle(Title[i]); 
      gvt[i]->GetXaxis()->SetTitle("Time [sec]"); 
      gvt[i]->GetXaxis()->CenterTitle(); 
      gvt[i]->GetYaxis()->SetTitle(Title[i]); 
      gvt[i]->GetYaxis()->CenterTitle(); 
   }

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(g[0],"SBS" ,"l");
   L->AddEntry(h[0],"LHRS","l");

   TCanvas *c1a = new TCanvas("c1a","BCM Check: LHRS and SBS",1200,800);
   c1a->Divide(2,2);

   TCanvas *c1b = new TCanvas("c1b","BCM Check: LHRS and SBS",1200,800);
   c1b->Divide(2,2);

   for(int i=0;i<N/2;i++){
      c1a->cd(i+1);
      if(logScale) gPad->SetLogy(true);
      g[i]->Draw();
      if(sameCanvas) h[i]->Draw("same");
      if(sameCanvas) L->Draw("same");
      c1a->Update();
      c1b->cd(i+1);
      if(logScale) gPad->SetLogy(true);
      g[i+4]->Draw();
      if(sameCanvas) h[i+4]->Draw("same");
      if(sameCanvas) L->Draw("same");
      c1b->Update();
   }

   TCanvas *c2a = new TCanvas("c2a","BCM Check: LHRS and SBS",1200,800);
   c2a->Divide(2,2);

   TCanvas *c2b = new TCanvas("c2b","BCM Check: LHRS and SBS",1200,800);
   c2b->Divide(2,2);

   for(int i=0;i<N/2;i++){
      c2a->cd(i+1);
      if(logScale) gPad->SetLogy(true);
      gvt[i]->Draw();
      if(sameCanvas) hvt[i]->Draw("same");
      if(sameCanvas) L->Draw("same");
      c2a->Update();
      c2b->cd(i+1);
      if(logScale) gPad->SetLogy(true);
      gvt[i+4]->Draw();
      if(sameCanvas) hvt[i+4]->Draw("same");
      if(sameCanvas) L->Draw("same");
      c2b->Update();
   }

   // test 60 Hz and 125 MHz clocks 
   // 60 Hz (last channel, "225". I think it's actually 223.) 
   TSsbs->Project("g60(100,0,100)"          ,"sbs.7_1.scalerRate","","");
   // 125 Hz (absolute channel 176, 208)  
   TSsbs->Project("g125a(100,124.8E+6,125.5E+6)","sbs.5_16.125MHz_CLK.rate","","");
   TSsbs->Project("g125b(100,124.8E+6,125.5E+6)","sbs.6_16.125MHz_CLK.rate","","");

   TH1F **gclk = new TH1F*[3]; 
   gclk[0] = (TH1F *)gDirectory->Get("g60"  ); 
   gclk[1] = (TH1F *)gDirectory->Get("g125a"); 
   gclk[2] = (TH1F *)gDirectory->Get("g125b"); 

   TCanvas *c3 = new TCanvas("c3","Additional clocks",1200,500);
   c3->Divide(3,1); 

   for(int i=0;i<3;i++){
      c3->cd(i+1);
      gclk[i]->Draw(); 
      c3->Update(); 
   }

   return 0;
}
//______________________________________________________________________________
int GetEvents(int run,std::string system,TTree *T){

   double clk1c=0,clk1r=0,clk2c=0,clk2r=0; // clk1 = 103.7 kHz, clk2 = 1024 Hz
   double u1r=-1,u1c=-1;
   double d1r=-1,d1c=-1;
   double d3r=-1,d3c=-1;
   double d10r=-1,d10c=-1;
   double dnewr=-1,dnewc=-1;
   double unewr=-1,unewc=-1;
   double unsr=-1,unsc=-1;
   
   std::cout << __LINE__ << std::endl;

   T->SetBranchAddress( Form("%s.bcm.u1.cnt"    ,system.c_str()) ,&u1c  );
   T->SetBranchAddress( Form("%s.bcm.u1.rate"   ,system.c_str()) ,&u1r  );
   T->SetBranchAddress( Form("%s.bcm.d1.cnt"    ,system.c_str()) ,&d1c  );
   T->SetBranchAddress( Form("%s.bcm.d1.rate"   ,system.c_str()) ,&d1r  );
   T->SetBranchAddress( Form("%s.bcm.d3.cnt"    ,system.c_str()) ,&d3c  );
   T->SetBranchAddress( Form("%s.bcm.d3.rate"   ,system.c_str()) ,&d3r  );
   T->SetBranchAddress( Form("%s.bcm.d10.cnt"   ,system.c_str()) ,&d10c );
   T->SetBranchAddress( Form("%s.bcm.d10.rate"  ,system.c_str()) ,&d10r );
   T->SetBranchAddress( Form("%s.bcm.dnew.cnt"  ,system.c_str()) ,&dnewc);
   T->SetBranchAddress( Form("%s.bcm.dnew.rate" ,system.c_str()) ,&dnewr);
   T->SetBranchAddress( Form("%s.bcm.unew.cnt"  ,system.c_str()) ,&unewc);
   T->SetBranchAddress( Form("%s.bcm.unew.rate" ,system.c_str()) ,&unewr);
   T->SetBranchAddress( Form("%s.bcm.unser.cnt" ,system.c_str()) ,&unsc );
   T->SetBranchAddress( Form("%s.bcm.unser.rate",system.c_str()) ,&unsr );

   std::cout << __LINE__ << std::endl;

   if(system.compare("Left")==0){
      T->SetBranchAddress("Left.slot02.T8.cnt"       ,&clk1c);
      T->SetBranchAddress("Left.slot02.RCLK1024.cnt" ,&clk2c);
      T->SetBranchAddress("Left.slot02.T8.rate"      ,&clk1r);
      T->SetBranchAddress("Left.slot02.RCLK1024.cnt" ,&clk2c);
      T->SetBranchAddress("Left.slot02.RCLK1024.rate",&clk2r);
   }
   // FIXME: What slots/channels are the clocks in for sbsvme29? 
 
   // char outpath[200];
   // sprintf(outpath,"rf_dump_%s_run-%d.csv",system.c_str(),run);

   // std::ofstream outfile;
   // outfile.open(outpath);
   // if(outfile.fail()){
   //    std::cout << "Cannot open the file: " << outpath << std::endl;
   // }

   // std::string header = "unewc,unewr,dnewc,dnewr,unsc,unsr,u1c,u1r,d1c,d1r,d3c,d3r,d10c,d10r,clk_103700_c,clk_103700_r,clk_1024_c,clk_1024_r";
   // outfile << header << std::endl;

   char msg[250];
   
   std::cout << __LINE__ << std::endl;
 
   int NEntries = T->GetEntries();
   for(int i=0;i<NEntries;i++){
      T->GetEntry(i);
      sprintf(msg,"%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf,%.0lf,%.3lf",
           	   unewc,unewr,dnewc,dnewr,unsc,unsr,u1c,u1r,d1c,d1r,d3c,d3r,d10c,d10r,clk1c,clk1r,clk2c,clk2r);
      // outfile << msg << std::endl;
      if(gIsDebug) std::cout << msg << std::endl;
   } 

   // outfile.close();
   // std::cout << "Data written to file: " << outpath << std::endl;
   
   return 0;
}

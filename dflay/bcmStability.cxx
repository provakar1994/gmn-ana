// Plot all BCM data vs run number  

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"

#include "./include/codaRun.h"
#include "./src/JSONManager.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

int bcmStability(const char *confPath){

   gStyle->SetOptStat(0);

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix  = jmgr->GetValueFromKey_str("ROOTfile-path"); 
   std::string runPath = jmgr->GetValueFromKey_str("run-path");  
   std::string cutPath = jmgr->GetValueFromKey_str("cut-path"); 
   std::string ccPath  = jmgr->GetValueFromKey_str("bcm-cc-path"); 
   std::string varType = jmgr->GetValueFromKey_str("var-type"); // rate, cnt, or current?  

   std::vector<codaRun_t> runList;  
   rc = bcm_util::LoadRuns(runPath.c_str(),runList);
   if(rc!=0) return 1; 

   BCMManager *mgr = new BCMManager("NONE",ccPath.c_str(),false);

   TString filePath;  
   const int NR = runList.size();  
   for(int i=0;i<NR;i++){ 
      filePath = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d.root",
                      prefix.c_str(),runList[i].runNumber,runList[i].stream,runList[i].segmentBegin,runList[i].segmentEnd);
      mgr->LoadFile(filePath,runList[i].runNumber);
   }

   // do stats by run number 
   std::vector<scalerData_t> rawData,data;
   mgr->GetVector_scaler("sbs",rawData);  

   // load cuts 
   std::vector<cut_t> cutList; 
   cut_util::LoadCuts(cutPath.c_str(),cutList);
   // apply cuts  
   cut_util::ApplyCuts(cutList,rawData,data);

   // sort the data by run number 
   std::sort(data.begin(),data.end(),compareScalerData_byRun); 

   const int N = 7;
   TString var[N]; 
   TString varRoot[N] = {"unser","u1","unew","d1","d3","d10","dnew"}; 
   for(int i=0;i<N;i++) var[i] = Form("%s.%s",varRoot[i].Data(),varType.c_str());  

   TGraphErrors **g  = new TGraphErrors*[N]; 

   std::vector<double> run;
   std::vector<double> mean_uns,stdev_uns;
   std::vector<double> mean_unew,stdev_unew;
   std::vector<double> mean_u1,stdev_u1;
   std::vector<double> mean_d1,stdev_d1;
   std::vector<double> mean_d3,stdev_d3;
   std::vector<double> mean_d10,stdev_d10;
   std::vector<double> mean_dnew,stdev_dnew;

   // get run stats for all variables  
   bcm_util::GetStats_byRun(var[0].Data(),data,run,mean_uns ,stdev_uns );
   run.clear();
   bcm_util::GetStats_byRun(var[1].Data(),data,run,mean_u1  ,stdev_u1  );
   run.clear();
   bcm_util::GetStats_byRun(var[2].Data(),data,run,mean_unew,stdev_unew);
   run.clear();
   bcm_util::GetStats_byRun(var[3].Data(),data,run,mean_d1  ,stdev_d1  );
   run.clear();
   bcm_util::GetStats_byRun(var[4].Data(),data,run,mean_d3  ,stdev_d3  );
   run.clear();
   bcm_util::GetStats_byRun(var[5].Data(),data,run,mean_d10 ,stdev_d10 );
   run.clear();
   bcm_util::GetStats_byRun(var[6].Data(),data,run,mean_dnew,stdev_dnew);

   // if we're missing runs in the initial run list, 
   // the analyzed run list is different (smaller) 
   int NNR = run.size();

   // print to screen
   for(int i=0;i<NNR;i++){
//       std::cout << Form("run %05d: unser = %.3lf ± %.3lf, u1 = %.3lf ± %.3lf, unew = %.3lf ± %.3lf, d1 = %.3lf ± %.3lf, d3 = %.3lf ± %.3lf, d10 = %.3lf ± %.3lf, dnew = %.3lf ± %.3lf",
// 	    (int)run[i],mean_uns[i],stdev_uns[i],mean_u1[i],stdev_u1[i],mean_unew[i],stdev_unew[i],mean_d1[i],stdev_d1[i],mean_d3[i],stdev_d3[i],mean_d10[i],stdev_d10[i],mean_dnew[i],stdev_dnew[i]) << std::endl;
      std::cout << Form("%05d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
	                (int)run[i],mean_uns[i],stdev_uns[i],mean_u1[i],stdev_u1[i],
                        mean_unew[i],stdev_unew[i],mean_d1[i],stdev_d1[i],mean_d3[i],stdev_d3[i],
                        mean_d10[i],stdev_d10[i],mean_dnew[i],stdev_dnew[i]) << std::endl;
   }

   // track unser change run to run
   double uns_thr = 2E+3; // 2 kHz => ~0.5 uA  
   double delta=0,err=0; 
   std::vector<int> runMarker; 
   std::vector<double> deltaUnser,deltaUnserErr; 
   for(int i=0;i<NNR;i++){
      if(i==0){
	 // first run in the list 
	 delta = 0;
	 err   = stdev_uns[i];
      }else{
	 delta = mean_uns[i] - mean_uns[i-1];
	 err   = TMath::Sqrt( stdev_uns[i]*stdev_uns[i] + stdev_uns[i-1]*stdev_uns[i-1] );
      }
      // if delta > threshold, this marks a new calibration period 
      if(delta>uns_thr){
	 // std::cout << "Calibration set for run " << run[i] << std::endl;
	 runMarker.push_back(run[i]); 
      } 
      deltaUnser.push_back(delta); 
      deltaUnserErr.push_back(err);
   }

   TGraphErrors *gDeltaUnser = graph_df::GetTGraphErrors(run,deltaUnser,deltaUnserErr); 
   graph_df::SetParameters(gDeltaUnser,20,kBlack); 

   int M=0;
   for(int i=0;i<N;i++){
      g[i] = bcm_util::GetTGraphErrors_byRun(var[i].Data(),data);
      graph_df::SetParameters(g[i],20,kBlack);
   }

   TGraphErrors *gUnserCurrent_byRun = bcm_util::GetTGraphErrors_byRun("unser.current",data);
   graph_df::SetParameters(gUnserCurrent_byRun,20,kBlue);

   TGraph *gUnserCurrent_byEvent = mgr->GetTGraph("sbs","event","unser.current");  
   graph_df::SetParameters(gUnserCurrent_byEvent,20,kBlue);

   // make a line at 0 uA 

   int runDelta = 50; 
   int runMin = run[0] - runDelta; 
   int runMax = run[NNR-1] + runDelta; 
   TLine *zero = new TLine(runMin,0,runMax,0); 
   zero->SetLineColor(kBlack); 

   TString Title,yAxisTitle,units;
   TString xAxisTitle = Form("Run Number");

   if(varType.compare("rate")==0){
      units = Form("Hz");
   }else if(varType.compare("current")==0){
      units = Form("#muA");
   }

   TCanvas *c1a = new TCanvas("c1a","BCM Check by Run",1200,800);
   c1a->Divide(2,2);

   TCanvas *c1b = new TCanvas("c1b","BCM Check by Run",1200,800);
   c1b->Divide(2,2); 

   for(int i=0;i<N/2;i++){
      c1a->cd(i+1);
      Title      = Form("%s"     ,var[i].Data());
      yAxisTitle = Form("%s [%s]",var[i].Data(),units.Data());
      g[i]->Draw("alp");
      graph_df::SetLabels(g[i],Title,xAxisTitle,yAxisTitle);
      g[i]->Draw("alp");
      c1a->Update();
      // next canvas 
      c1b->cd(i+1);
      Title      = Form("%s"     ,var[i+3].Data());
      yAxisTitle = Form("%s [%s]",var[i+3].Data(),units.Data());
      g[i+3]->Draw("alp");
      graph_df::SetLabels(g[i+3],Title,xAxisTitle,yAxisTitle);
      g[i+3]->Draw("alp");
      c1b->Update();
   }

   // last one 
   Title      = Form("%s"     ,var[6].Data());
   yAxisTitle = Form("%s [%s]",var[6].Data(),units.Data());
   c1b->cd(4);
   g[6]->Draw("alp");
   graph_df::SetLabels(g[6],Title,xAxisTitle,yAxisTitle);
   g[6]->Draw("alp");
   c1b->Update();

   xAxisTitle = Form("Unser Current [#muA]"); 

   TCanvas *c2 = new TCanvas("c2","Unser Current",1200,800);
   c2->Divide(1,2);
 
   c2->cd(1);
   gUnserCurrent_byEvent->Draw("ap");
   graph_df::SetLabels(gUnserCurrent_byEvent,"Unser Current","Event","Unser Current [#muA]");
   gUnserCurrent_byEvent->Draw("ap");
   c2->Update(); 

   c2->cd(2);
   gUnserCurrent_byRun->Draw("ap");
   graph_df::SetLabels(gUnserCurrent_byRun,"Unser Current","Run Number","Unser Current [#muA]");
   gUnserCurrent_byRun->GetXaxis()->SetLimits(runMin,runMax);
   gUnserCurrent_byRun->Draw("ap");
   zero->Draw("same"); 
   c2->Update(); 
   
   TCanvas *c3 = new TCanvas("c3","Change in Unser",1200,600); 
   c3->Divide(1,2);

   c3->cd(1);
   g[0]->Draw("alp");
   graph_df::SetLabels(g[0],"Unser Rate","Run Number","Unser Rate [Hz]");  
   g[0]->Draw("alp");
   c3->Update();

   c3->cd(2);
   gDeltaUnser->Draw("alp");
   graph_df::SetLabels(gDeltaUnser,"Change in Unser","Run Number","#Delta Unser Rate [Hz]");  
   gDeltaUnser->Draw("alp");
   c3->Update();

   return 0;
}

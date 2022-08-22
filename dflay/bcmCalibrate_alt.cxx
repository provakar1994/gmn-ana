// Plot all BCM data vs Unser current 
// This is the "in-run" method  

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"

#include "./include/codaRun.h"
#include "./src/Graph.cxx"
#include "./src/BCMManager.cxx"
#include "./src/CSVManager.cxx"
#include "./src/JSONManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

double fitFunc(double *x,double *p); 

int bcmCalibrate_alt(const char *confPath){

   // settings 
   bool logScale   = false;

   const int N = 7;
   TString varRoot[N] = {"unser","u1","unew","d1","d3","d10","dnew"};
   TString var[N]; 
   for(int i=0;i<N;i++) var[i] = Form("%s.rate.ps",varRoot[i].Data());  

   gStyle->SetOptStat(0);

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix    = jmgr->GetValueFromKey_str("ROOTfile-path"); 
   std::string runPath   = jmgr->GetValueFromKey_str("run-path");  
   std::string cutPath   = jmgr->GetValueFromKey_str("cut-path");  
   std::string bcmCCPath = jmgr->GetValueFromKey_str("bcm-cc-path");  
   std::string outdir    = jmgr->GetValueFromKey_str("output-dir");  

   std::vector<codaRun_t> runList;  
   rc = util_df::LoadRunList(runPath.c_str(),prefix.c_str(),runList);
   if(rc!=0) return 1; 

   BCMManager *mgr = new BCMManager("NONE",bcmCCPath.c_str(),false);

   rc = util_df::LoadBCMData(runList,mgr);
   if(rc!=0) return 1; 

   // do stats by run number 
   std::vector<scalerData_t> rawData,data;
   mgr->GetVector_scaler("sbs",rawData);  

   // load cuts (loose cut on the Unser)  
   std::vector<cut_t> cutList; 
   cut_util::LoadCuts(cutPath.c_str(),cutList);
   // apply cuts  
   cut_util::ApplyCuts(cutList,rawData,data);

   // sort the data by run number 
   std::sort(data.begin(),data.end(),compareScalerData_byRun); 
   
   std::vector<int> run;
   mgr->GetRunList(run);
   const int NNR = run.size();
   if(NNR==0){
      std::cout << "ERROR! Run list size = " << NNR << std::endl;
      return 1;
   }

   // we load the requested beam current into the codaRun::info data field as a check
   // so we will populate this vector 
   // double argC=0; 
   std::vector<double> reqCurrent,reqCurrentErr; 
   // for(int i=0;i<NNR;i++){
   //    for(int j=0;j<NR;j++){
   //       if(run[i]==runList[j].runNumber){
   //          argC = std::atof(runList[j].info.c_str());
   //          // reqCurrent.push_back(argC); 
   //       }
   //    }
   // } 

   std::vector<epicsData_t> eRAW_DATA,eDATA;
   std::vector<scalerData_t> RAW_DATA,DATA;
   std::vector<cut_t> runCuts,erunCuts; 

   char inpath_cuts[200]; 

   std::vector<double> mean_ibc,stdev_ibc; 
   std::vector<double> mean_hac,stdev_hac; 

   std::vector<double> rr,mean_uns,stdev_uns;
   std::vector<double> mean_unsi,stdev_unsi;
   std::vector<double> diff,diffErr; // difference between unser current and requested current 

   std::vector<double> RR,UR,URE,UI,UIE; 
   std::vector<double> U1,U1E,UN,UNE;
   std::vector<double> D1,D1E,DN,DNE;
   std::vector<double> D3,D3E,D10,D10E;
   std::vector<double> hacCurrent,hacCurrentErr;

   std::vector<double> mean_u1,stdev_u1;
   std::vector<double> mean_unew,stdev_unew;
   std::vector<double> mean_d1,stdev_d1;
   std::vector<double> mean_d3,stdev_d3;
   std::vector<double> mean_d10,stdev_d10;
   std::vector<double> mean_dnew,stdev_dnew;

   int NEV=0,N_CUTS=0; 
   double argDiff=0,argDiffErr=0;
   for(int i=0;i<NNR;i++){
      std::cout << Form("Processing run %d",(int)run[i]) << std::endl;
      // grab a run and apply a list of cuts 
      mgr->GetVector_scaler("sbs",run[i],RAW_DATA);
      mgr->GetVector_epics(run[i],eRAW_DATA);
      // load and apply cuts
      sprintf(inpath_cuts,"./input/cut-list/bcm-calib/run-%d.json",run[i]);  
      cut_util::LoadCuts_json(inpath_cuts,runCuts);
      cut_util::LoadCuts_epics_json(inpath_cuts,erunCuts);
      N_CUTS = runCuts.size();
      if(N_CUTS>0){
	 cut_util::ApplyCuts_alt(runCuts,RAW_DATA,DATA);
      }else{
	 // copy all data, no cuts
	 NEV = RAW_DATA.size();
	 for(int j=0;j<NEV;j++) DATA.push_back(RAW_DATA[j]); 
      } 
      // epics data 
      N_CUTS = erunCuts.size();
      if(N_CUTS>0){
	 cut_util::ApplyCuts_alt(erunCuts,eRAW_DATA,eDATA);
      }else{
	 // copy all data, no cuts
	 NEV = eRAW_DATA.size();
	 for(int j=0;j<NEV;j++) eDATA.push_back(eRAW_DATA[j]);
      }
      // do some stats 
      bcm_util::GetStats_byRun("unser.current",DATA,rr,mean_unsi,stdev_unsi);
      rr.clear();
      bcm_util::GetStats_byRun(var[0].Data(),DATA,rr,mean_uns,stdev_uns);
      // get run stats for all variablData
      bcm_util::GetStats_byRun(var[1].Data(),DATA,rr,mean_u1  ,stdev_u1  );
      rr.clear();
      bcm_util::GetStats_byRun(var[2].Data(),DATA,rr,mean_unew,stdev_unew);
      rr.clear();
      bcm_util::GetStats_byRun(var[3].Data(),DATA,rr,mean_d1  ,stdev_d1  );
      rr.clear();
      bcm_util::GetStats_byRun(var[4].Data(),DATA,rr,mean_d3  ,stdev_d3  );
      rr.clear();
      bcm_util::GetStats_byRun(var[5].Data(),DATA,rr,mean_d10 ,stdev_d10 );
      rr.clear();
      bcm_util::GetStats_byRun(var[6].Data(),DATA,rr,mean_dnew,stdev_dnew);
      // epics current 
      bcm_util::GetStats_byRun_epics("IBC1H04CRCUR2"  ,eDATA,rr,mean_ibc,stdev_ibc); 
      rr.clear();
      bcm_util::GetStats_byRun_epics("hac_bcm_average",eDATA,rr,mean_hac,stdev_hac);
      // print to screen
      std::cout << Form("run %d: EPICS IBC1H04CRCUR2 = %.3lf ± %.3lf, EPICS hac_bcm_average = %.3lf ± %.3lf, unser current = %.3lf ± %.3lf",
                        (int)rr[0],mean_ibc[0],stdev_ibc[0],mean_hac[0],stdev_hac[0],mean_unsi[0],stdev_unsi[0]) << std::endl;
      // save for diagnostic plot
      reqCurrent.push_back(mean_ibc[0]);  
      reqCurrentErr.push_back(stdev_ibc[0]);  
      hacCurrent.push_back(mean_hac[0]);  
      hacCurrentErr.push_back(stdev_hac[0]);  
      argDiff    = (mean_ibc[0] - mean_unsi[0]);
      argDiffErr = TMath::Sqrt( stdev_ibc[0]*stdev_ibc[0] + stdev_unsi[0]*stdev_unsi[0]);  
      diff.push_back(argDiff);
      diffErr.push_back(argDiffErr);  
      RR.push_back(rr[0]);
      UR.push_back(mean_uns[0]); 
      URE.push_back(stdev_uns[0]);  
      UI.push_back(mean_unsi[0]); 
      UIE.push_back(stdev_unsi[0]);
      U1.push_back(mean_u1[0]); 
      U1E.push_back(stdev_u1[0]);  
      UN.push_back(mean_unew[0]); 
      UNE.push_back(stdev_unew[0]);  
      D1.push_back(mean_d1[0]); 
      D1E.push_back(stdev_d1[0]);  
      D3.push_back(mean_d3[0]); 
      D3E.push_back(stdev_d3[0]);  
      D10.push_back(mean_d10[0]); 
      D10E.push_back(stdev_d10[0]);  
      DN.push_back(mean_dnew[0]); 
      DNE.push_back(stdev_dnew[0]);  
      // clear  
      rr.clear();
      mean_ibc.clear();
      stdev_ibc.clear();
      mean_hac.clear();
      stdev_hac.clear();
      mean_uns.clear();
      mean_unsi.clear();
      stdev_uns.clear(); 
      stdev_unsi.clear();
      mean_u1.clear();
      stdev_u1.clear();
      mean_unew.clear();
      stdev_unew.clear();
      mean_d1.clear();
      stdev_d1.clear(); 
      mean_d3.clear();
      stdev_d3.clear();
      mean_d10.clear();
      stdev_d10.clear();
      mean_dnew.clear();
      stdev_dnew.clear();
      // set up for next run
      RAW_DATA.clear();
      DATA.clear();
      runCuts.clear(); 
      eRAW_DATA.clear();
      eDATA.clear();
      erunCuts.clear(); 
   }
   std::cout << "----------------------------------" << std::endl;

   TGraphErrors *gReqCurrentRun   = graph_df::GetTGraphErrors(RR,reqCurrent,reqCurrentErr);
   TGraphErrors *gHACCurrentRun   = graph_df::GetTGraphErrors(RR,hacCurrent,hacCurrentErr);
   TGraphErrors *gUnserCurrentRun = graph_df::GetTGraphErrors(RR,UI,UIE);
   TGraphErrors *gUnserRateRun    = graph_df::GetTGraphErrors(RR,UR,URE);
   TGraphErrors *gUnserDiffRun    = graph_df::GetTGraphErrors(RR,diff,diffErr); 
   graph_df::SetParameters(gReqCurrentRun  ,21,kBlack); 
   graph_df::SetParameters(gHACCurrentRun  ,22,kViolet); 
   graph_df::SetParameters(gUnserCurrentRun,20,kRed); 
   graph_df::SetParameters(gUnserRateRun   ,20,kRed);
   graph_df::SetParameters(gUnserDiffRun   ,20,kBlue);

   // make a line at 0 uA 
   int runDelta = 50; 
   int runMin = run[0] - runDelta; 
   int runMax = run[NNR-1] + runDelta; 
   TLine *zero = new TLine(runMin,0,runMax,0); 
   zero->SetLineColor(kBlack); 

   TMultiGraph *mgu = new TMultiGraph();
   mgu->Add(gReqCurrentRun  ,"p");
   mgu->Add(gHACCurrentRun  ,"p");
   mgu->Add(gUnserCurrentRun,"p");

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);
   L->AddEntry(gReqCurrentRun  ,"IBC1H04CRCUR2"  ,"p"); 
   L->AddEntry(gHACCurrentRun  ,"hac_bcm_average","p"); 
   L->AddEntry(gUnserCurrentRun,"Unser Current"  ,"p"); 

   TCanvas *c4 = new TCanvas("c4","Unser Data",1200,800); 
   c4->Divide(1,2);

   c4->cd(1); 
   mgu->Draw("a"); 
   graph_df::SetLabels(mgu,"Beam Current","Run Number","Beam Current [#muA]");  
   graph_df::SetLabelSizes(mgu,0.05,0.06);  
   mgu->GetXaxis()->SetLimits(runMin,runMax);  
   mgu->GetYaxis()->SetLimits(-2,20);  
   mgu->Draw("a");
   L->Draw("same"); 
   c4->Update();

   c4->cd(2); 
   gUnserDiffRun->Draw("ap"); 
   graph_df::SetLabels(gUnserDiffRun,"Difference (IBC1H04CRCUR2 - Unser)","Run Number","Diff (EPICS - Unser) [#muA]"); 
   graph_df::SetLabelSizes(gUnserDiffRun,0.05,0.06);  
   gUnserDiffRun->GetXaxis()->SetLimits(runMin,runMax);  
   gUnserDiffRun->GetYaxis()->SetRangeUser(-2,2);  
   gUnserDiffRun->Draw("ap");
   zero->Draw("same"); 
   c4->Update();

   // build the scaler rate vs Unser current plots
   TGraphErrors **g  = new TGraphErrors*[N];
   g[0] = graph_df::GetTGraphErrors(UI,UIE,UR,URE); 
   g[1] = graph_df::GetTGraphErrors(UI,UIE,U1,U1E); 
   g[2] = graph_df::GetTGraphErrors(UI,UIE,UN,UNE); 
   g[3] = graph_df::GetTGraphErrors(UI,UIE,D1,D1E); 
   g[4] = graph_df::GetTGraphErrors(UI,UIE,D3,D3E); 
   g[5] = graph_df::GetTGraphErrors(UI,UIE,D10,D10E); 
   g[6] = graph_df::GetTGraphErrors(UI,UIE,DN,DNE); 
 
   for(int i=0;i<N;i++){
      graph_df::SetParameters(g[i],20,kBlack);
   }

   // plot the scaler rates vs Unser current 
   TString Title,yAxisTitle;
   TString xAxisTitle = Form("Unser Current [#muA]");

   TCanvas *c1a = new TCanvas("c1a","BCMs vs Unser Current",1200,800);
   c1a->Divide(2,2);

   TCanvas *c1b = new TCanvas("c1b","BCMs vs Unser Current",1200,800);
   c1b->Divide(2,2); 

   double par[N][2],parErr[N][2]; 

   const int npar=2;
   // double min=-1,max=20; 
   // fit range     unser, u1, unew, d1, d3, d10, dnew
   double min[N] = {-1, 3,-1, 3, 3, 1,-1};
   double max[N] = {16,16,16,16,16,16,16};
   TString fitName[N]; 
   double argPar=0,argParErr=0;
   TF1 **myFit = new TF1*[N]; 
   for(int i=0;i<N;i++){
      fitName[i] = Form("myFit_%d",i);
      myFit[i] = new TF1(fitName[i],fitFunc,min[i],max[i],npar);
      myFit[i]->SetParameter(0,0); 
      myFit[i]->SetParameter(1,0); 
      myFit[i]->SetLineColor(kRed);
   }

   for(int i=0;i<N/2;i++){
      c1a->cd(i+1);
      gStyle->SetOptFit(111); 
      Title      = Form("%s"     ,var[i].Data());
      yAxisTitle = Form("%s [Hz]",var[i].Data());
      g[i]->Draw("ap");
      graph_df::SetLabels(g[i],Title,xAxisTitle,yAxisTitle);
      g[i]->Draw("ap");
      g[i]->Fit(fitName[i],"QR"); 
      for(int j=0;j<npar;j++){
	 argPar    = myFit[i]->GetParameter(j);
	 argParErr = myFit[i]->GetParError(j);
         par[i][j] = argPar; 
         parErr[i][j] = argParErr; 
      }
      c1a->Update();
      // next canvas 
      c1b->cd(i+1);
      gStyle->SetOptFit(111); 
      Title      = Form("%s"     ,var[i+3].Data());
      yAxisTitle = Form("%s [Hz]",var[i+3].Data());
      g[i+3]->Draw("ap");
      graph_df::SetLabels(g[i+3],Title,xAxisTitle,yAxisTitle);
      g[i+3]->Draw("ap");
      g[i+3]->Fit(fitName[i+3],"QR"); 
      for(int j=0;j<npar;j++){
	 argPar    = myFit[i+3]->GetParameter(j);
	 argParErr = myFit[i+3]->GetParError(j);
         par[i+3][j] = argPar; 
         parErr[i+3][j] = argParErr; 
      }
      c1b->Update();
   }

   // last one 
   Title      = Form("%s"     ,var[6].Data());
   yAxisTitle = Form("%s [Hz]",var[6].Data());
   c1b->cd(4);
   gStyle->SetOptFit(111); 
   g[6]->Draw("ap");
   graph_df::SetLabels(g[6],Title,xAxisTitle,yAxisTitle);
   g[6]->Draw("ap");
   g[6]->Fit(fitName[6],"QR"); 
   for(int j=0;j<npar;j++){
      argPar    = myFit[6]->GetParameter(j);
      argParErr = myFit[6]->GetParError(j);
      par[6][j] = argPar; 
      parErr[6][j] = argParErr; 
   }
   c1b->Update();
  
   // gather results into calibCoeff struct and print to file  
   calibCoeff_t apt; 
   std::vector<calibCoeff_t> cc; 

   // create the output directory 
   util_df::MakeDirectory(outdir.c_str());

   // write to individual files
   char outpath[200]; 
   for(int i=0;i<N;i++){
      apt.dev       = varRoot[i].Data(); 
      apt.offset    = par[i][0]; 
      apt.offsetErr = parErr[i][0];  
      apt.slope     = par[i][1]; 
      apt.slopeErr  = parErr[i][1]; 
      apt.runMin    = 0; 
      apt.runMax    = 99999;
      std::cout << Form("%s: offset = %.2lf ± %.2lf Hz, slope = %.2lf ± %.2lf",
                        apt.dev.c_str(),apt.offset,apt.offsetErr,apt.slope,apt.slopeErr) << std::endl; 
      // file name 
      sprintf(outpath,"%s/%s-calib-results.csv",outdir._str(),apt.dev.c_str());
      cc.push_back(apt);
      bcm_util::WriteToFile_cc(outpath,cc);
      // set up for next device 
      cc.clear(); 
   }
 
   return 0;
}
//______________________________________________________________________________
double fitFunc(double *x,double *p){
   double f = p[0] + p[1]*x[0]; 
   return f;
}

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
#include "./src/CSVManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/JSONManager.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

TString GetTitle(std::vector<int> run); 

int GetStats(std::vector<std::string> var,std::vector<scalerData_t> data,
             std::vector<double> &run,std::vector<std::vector<double>> &mean,std::vector<std::vector<double>> &stdev);

int GetCharge(std::vector<std::string> var,BCMManager *mgr,
              std::vector<double> &run,std::vector<std::vector<double>> &mean,std::vector<std::vector<double>> &stdev); 

int GetCharge(std::string var,std::vector<scalerData_t> allData,std::vector<double> run,
              std::vector<double> &Q,std::vector<double> &dQ,std::vector<double> &Time,std::vector<double> &dTime); 

double getChargeError(double Q,double I, double dI,double t,double dt); 

int bcmCheckRun(const char *confPath){

   // settings 
   bool logScale = false;

   gStyle->SetOptStat(0);

   int rc=0;

   // read input configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string prefix       = jmgr->GetValueFromKey_str("ROOTfile-path"); 
   std::string runPath      = jmgr->GetValueFromKey_str("run-path");  
   std::string cutPath      = jmgr->GetValueFromKey_str("cut-path"); 
   std::string bcmCalibPath = jmgr->GetValueFromKey_str("bcm-cc-path");   
   
   BCMManager *mgr = new BCMManager("NONE",bcmCalibPath.c_str(),false);

   std::vector<codaRun_t> runList;  
   rc = util_df::LoadRunList(runPath.c_str(),prefix.c_str(),runList);
   if(rc!=0) return 1; 

   util_df::LoadBCMData(runList,mgr); 
   if(rc!=0) return 1; 

   // do stats by run number 
   std::vector<scalerData_t> rawData,data,RAW_DATA,DATA;
   mgr->GetVector_scaler("sbs",rawData);  

   std::vector<epicsData_t> erawData,edata,eRAW_DATA,eDATA; 
   mgr->GetVector_epics(edata);  

   // load cuts 
   std::vector<cut_t> cutList; 
   cut_util::LoadCuts(cutPath.c_str(),cutList);
   // apply cuts  
   cut_util::ApplyCuts(cutList,rawData,data);

   // sort the data by run number 
   std::sort(data.begin(),data.end(),compareScalerData_byRun); 
   std::sort(edata.begin(),edata.end(),compareEPICSData_byRun); 

   std::vector<int> rr; 
   mgr->GetRunList(rr); 
   const int NNR = rr.size();

   char inpath_cuts[200]; 
   std::vector<cut_t> runCuts,erunCuts; 

   // load and apply cuts
   int NEV=0,N_CUTS=0; 
   for(int i=0;i<NNR;i++){
      // grab a run and apply a list of cuts 
      mgr->GetVector_scaler("sbs",rr[i],RAW_DATA);
      mgr->GetVector_epics(rr[i],eRAW_DATA);
      // load and apply cuts 
      sprintf(inpath_cuts,"./input/cut-list/bcm-calib/run-%d.json",rr[i]);
      cut_util::LoadCuts_json(inpath_cuts,runCuts);
      cut_util::LoadCuts_epics_json(inpath_cuts,erunCuts);
      // scaler data 
      N_CUTS = runCuts.size();
      if(N_CUTS>0){
	 cut_util::ApplyCuts_alt(runCuts ,RAW_DATA ,DATA);
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
      // reset 
      runCuts.clear();
      erunCuts.clear();
      RAW_DATA.clear();  
      eRAW_DATA.clear();  
   }

   // NEV = DATA.size(); 
   // for(int i=0;i<NEV;i++) DATA[i].Print("rate"); 
   // return 0; 

   const int N = 7; 
   TString var[N]  = {"unser.rate"   ,"u1.rate"   ,"unew.rate"   ,"d1.rate"   ,"d3.rate"   ,"d10.rate"   ,"dnew.rate"};
   TString varC[N] = {"unser.current","u1.current","unew.current","d1.current","d3.current","d10.current","dnew.current"};

   // get the charge 
   // std::vector<double> RR;
   // std::vector<std::vector<double>> SS,MM; 
   // std::vector<std::string> VAR; 
   // for(int i=0;i<N;i++) VAR.push_back(varC[i].Data()); 
   // GetCharge(VAR,mgr,RR,MM,SS);

   std::vector<double> run,mean,stdev,q,dq,t,dt;
   std::vector<std::vector<double>> mm,ss,Q,dQ;     // no cuts
   std::vector<std::vector<double>> mmc,ssc,Qc,dQc; // with cuts 
   std::vector<std::vector<double>> mtc,stc;        // time of run 

   bcm_t dataPt; 
   std::vector<bcm_t> bcmData; 

   // get run stats for all variables
   int M=0;
   for(int i=0;i<N;i++){
      bcm_util::GetStats_byRun(varC[i].Data(),data,run,mean,stdev);
      GetCharge(varC[i].Data(),data,run,q,dq,t,dt); 
      // store results
      mm.push_back(mean); 
      ss.push_back(stdev);
      Q.push_back(q);  
      dQ.push_back(dq);  
      // now with cuts 
      run.clear();
      mean.clear();
      stdev.clear();
      q.clear();
      dq.clear();
      t.clear();
      dt.clear();
      bcm_util::GetStats_byRun(varC[i].Data(),DATA,run,mean,stdev);
      GetCharge(varC[i].Data(),DATA,run,q,dq,t,dt); 
      // store results
      mmc.push_back(mean); 
      ssc.push_back(stdev); 
      Qc.push_back(q);  
      dQc.push_back(dq);
      mtc.push_back(t); 
      stc.push_back(dt); 
      // dataPt.dev        = varC[i].Data();
      // dataPt.current    = mean; 
      // dataPt.currentErr = stdev; 
      // dataPt.charge     = q; 
      // dataPt.chargeErr  = dq; 
      // bcmData.push_back(dataPt);  
      // set up for next variable 
      run.clear();
      mean.clear();
      stdev.clear();
      q.clear();
      dq.clear();
      t.clear();
      dt.clear();
   } 

   // print results to screen
   for(int i=0;i<N;i++){
      std::cout << Form("%s: ",varC[i].Data()) << std::endl;
      M = mm[i].size(); // run dimension 
      for(int j=0;j<M;j++){
         std::cout << Form("   run %d: no cuts I = = %.3lf ± %.3lf, Q = %.5lf C, with cuts I = %.3lf ± %.3lf, Q = %.5lf C",
                           rr[j],mm[i][j],ss[i][j],Q[i][j],mmc[i][j],ssc[i][j],Qc[i][j]) << std::endl;
      }
   } 
  
   char msg[200]; 
 
   // print results to screen
   // beam current 
   M = mm[0].size(); 
   for(int i=0;i<M;i++){  // run dimension
      std::cout << "run " << rr[i] << ","; 
      for(int j=0;j<N;j++){ // BCM variable dimension
         if(j==0){
	    sprintf(msg,"%.3lf,%.3lf",mmc[j][i],ssc[j][i]); 
         }else{
	    sprintf(msg,"%s,%.3lf,%.3lf",msg,mmc[j][i],ssc[j][i]); 
         } 
      }
      std::cout << msg << std::endl;
   }

   // beam charge
   for(int i=0;i<M;i++){  // run dimension
      std::cout << "run " << rr[i] << ","; 
      for(int j=0;j<N;j++){ // BCM variable dimension
         dQc[j][i] = getChargeError(Qc[j][i],mmc[j][i],ssc[j][i],mtc[j][i],stc[j][i]);
         if(j==0){
	    sprintf(msg,"%.5lf,%.5lf",Qc[j][i],dQc[j][i]); 
         }else{
	    sprintf(msg,"%s,%.5lf,%.5lf",msg,Qc[j][i],dQc[j][i]); 
         } 
      }
      std::cout << msg << std::endl;
   }

   // get BCM plots vs event number
   std::string xVar = "time"; 
   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   int color[N] = {kBlack,kRed,kBlue,kGreen+2,kMagenta,kCyan+2,kOrange}; 
   TMultiGraph **mg  = new TMultiGraph*[N]; 
   TMultiGraph *mgc = new TMultiGraph(); 
   TGraph **gb  = new TGraph*[N]; 
   TGraph **ga  = new TGraph*[N]; 
   TGraph **gac = new TGraph*[N]; 
   for(int i=0;i<N;i++){
      mg[i]  = new TMultiGraph();
      gb[i]  = bcm_util::GetTGraph(xVar.c_str(),var[i].Data(),data); 
      ga[i]  = bcm_util::GetTGraph(xVar.c_str(),var[i].Data(),DATA); 
      gac[i] = bcm_util::GetTGraph(xVar.c_str(),varC[i].Data(),DATA); 
      graph_df::SetParameters(gb[i],21,kBlack);
      graph_df::SetParameters(ga[i],20,kRed);
      graph_df::SetParameters(gac[i],20,color[i]);
      mg[i]->Add(gb[i],"p"); 
      mg[i]->Add(ga[i],"p");
      mgc->Add(gac[i],"p"); 
      L->AddEntry(gac[i],varC[i].Data(),"p"); 
   }

   TString Title,yAxisTitle;
   // TString xAxisTitle = Form("Run Event Number");
   TString xAxisTitle = Form("Time");

   TString canvasTitle = GetTitle(rr); 
  
   TCanvas *c1a = new TCanvas("c1a",canvasTitle,1200,800);
   c1a->Divide(2,2);

   TCanvas *c1b = new TCanvas("c1b",canvasTitle,1200,800);
   c1b->Divide(2,2); 

   for(int i=0;i<N/2;i++){
      c1a->cd(i+1);
      Title      = Form("%s"     ,var[i].Data());
      yAxisTitle = Form("%s [Hz]",var[i].Data());
      mg[i]->Draw("alp");
      graph_df::SetLabels(mg[i],Title,xAxisTitle,yAxisTitle);
      graph_df::UseTimeDisplay(mg[i]); 
      mg[i]->Draw("alp");
      c1a->Update();
      // next canvas 
      c1b->cd(i+1);
      Title      = Form("%s"     ,var[i+3].Data());
      yAxisTitle = Form("%s [Hz]",var[i+3].Data());
      mg[i+3]->Draw("alp");
      graph_df::SetLabels(mg[i+3],Title,xAxisTitle,yAxisTitle);
      graph_df::UseTimeDisplay(mg[i+3]); 
      mg[i+3]->Draw("alp");
      c1b->Update();
   }

   // last one 
   Title      = Form("%s"     ,var[6].Data());
   yAxisTitle = Form("%s [Hz]",var[6].Data());
   c1b->cd(4);
   mg[6]->Draw("alp");
   graph_df::SetLabels(mg[6],Title,xAxisTitle,yAxisTitle);
   graph_df::UseTimeDisplay(mg[6]); 
   mg[6]->Draw("alp");
   c1b->Update();

   // make a plot of the beam current from EPICS
   const int NE = 2; 
   TString evar[NE] = {"IBC1H04CRCUR2","hac_bcm_average"};
   int ecolor[NE] = {kBlack,kRed};
   int eMarker[NE] = {20,21};

   TMultiGraph *mge = new TMultiGraph();

   TGraph **geb = new TGraph*[NE]; 
   TGraph **gea = new TGraph*[NE]; 
   
   TLegend *Le = new TLegend(0.6,0.6,0.8,0.8);

   std::vector<double> er,eCur,deCur; 
   std::vector<double> eCUR,deCUR; 

   for(int i=0;i<NE;i++){
      geb[i] = bcm_util::GetTGraph(xVar.c_str(),evar[i].Data(),edata);
      gea[i] = bcm_util::GetTGraph(xVar.c_str(),evar[i].Data(),eDATA);
      graph_df::SetParameters(geb[i],eMarker[i],kBlack);   
      graph_df::SetParameters(gea[i],eMarker[i],kRed);   
      mge->Add(geb[i],"p"); 
      mge->Add(gea[i],"p"); 
      Le->AddEntry(geb[i],evar[i],"p"); 
      // get run stats 
      bcm_util::GetStats_byRun_epics(evar[i].Data(),edata,er,eCur,deCur); 
      er.clear();
      bcm_util::GetStats_byRun_epics(evar[i].Data(),eDATA,er,eCUR,deCUR);
      std::cout << Form("EPICS %s:",evar[i].Data()) << std::endl;
      M = er.size();
      for(int j=0;j<M;j++){
	 std::cout << Form("   run %d: no cuts = %.3lf ± %.3lf, with cuts = %.3lf ± %.3lf",
                           (int)er[j],eCur[j],deCur[j],eCUR[j],deCUR[j]) << std::endl;   
      }
      // set up for next variable
      er.clear();
      eCur.clear();
      deCur.clear(); 
      eCUR.clear();
      deCUR.clear(); 
   }

   Title = GetTitle(rr); 

   TCanvas *c3 = new TCanvas("c3","Unser and BCM Current",1200,800);
   c3->Divide(1,2);

   c3->cd(1);
   mgc->Draw("a"); 
   graph_df::SetLabels(mgc,Title,xAxisTitle,"Beam Current [#muA]");
   graph_df::SetLabelSizes(mgc,0.05,0.06);  
   graph_df::UseTimeDisplay(mgc);  
   mgc->Draw("a");
   L->Draw("same"); 
   c3->Update(); 

   c3->cd(2);
   mge->Draw("a");
   graph_df::SetLabels(mge,"EPICS Current",xAxisTitle,"Beam Current [#muA]"); 
   graph_df::SetLabelSizes(mge,0.05,0.06); 
   // graph_df::UseTimeDisplay(mge);  
   mge->Draw("a");
   Le->Draw("same"); 
   c3->Update();

   delete mgr;
   delete jmgr; 

   return 0;
}
//______________________________________________________________________________
double getChargeError(double Q,double I, double dI,double t,double dt){
   double T1=0,T2=0;
   if(I!=0) T1 = TMath::Power(dI/I,2.); 
   if(t!=0) T2 = TMath::Power(dt/t,2.); 
   double dQ = Q*TMath::Sqrt( T1 + T2 );
   return dQ; 
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
//______________________________________________________________________________
int GetStats(std::vector<std::string> var,std::vector<scalerData_t> data,
             std::vector<double> &run,std::vector<std::vector<double>> &mean,std::vector<std::vector<double>> &stdev){
   // get mean and stdev for all BCMs for all runs  
   std::vector<double> mm,ss; 
   const int NV = var.size();
   for(int i=0;i<NV;i++){
      bcm_util::GetStats_byRun(var[i].c_str(),data,run,mm,ss);
      mean.push_back(mm);
      stdev.push_back(ss);
      // set up for next BCM variable
      mm.clear();
      ss.clear(); 
      if(i!=NV-1) run.clear(); // delete until last value  
   }
   return 0;
}
//______________________________________________________________________________
int GetCharge(std::string var,std::vector<scalerData_t> allData,std::vector<double> run,
              std::vector<double> &Q,std::vector<double> &dQ,std::vector<double> &Time,std::vector<double> &dTime){
   // get the charge associated with the run

   int M=0;
   const int NT = allData.size(); 
   const int NR = run.size();
   double MICROAMPS = 1E-6;  
   double timeStep=0,timeStep103=0,deltaTime=0,chargeSum=0,chargeSum103=0;
   std::vector<scalerData_t> runData;  
   for(int i=0;i<NR;i++){
      // get the run data
      for(int j=0;j<NT;j++){ 
	 if(allData[j].runNumber==(int)run[i]){
	    runData.push_back(allData[j]); 
         } 
      }
      // now loop over the run data
      M = runData.size();
      deltaTime = runData[M-1].time - runData[0].time;
      // accumulate charge over the whole run for each variable  
      for(int j=0;j<M;j++){
         if(j==0){
	    timeStep    = runData[1].time - runData[0].time;  
	    timeStep103 = runData[1].time - runData[0].time;  
	 }else{
	    timeStep    = runData[j].time - runData[j-1].time;
	    timeStep103 = runData[j].time - runData[j-1].time;
	 }
	 // if(j<10){
	 //    std::cout << Form("event %d, timeStep(RF time) = %.3lf sec, timeStep(103kHz) = %.3lf",
	 //                      j,timeStep,timeStep103) << std::endl;
	 // }
	 chargeSum    += timeStep*runData[j].getValue(var)*MICROAMPS;    // convert to Amps 
	 chargeSum103 += timeStep103*runData[j].getValue(var)*MICROAMPS; // convert to Amps 
      }
      // print to screen 
      std::cout << Form("run %d:",(int)run[i]) << std::endl;
      std::cout << Form("   %15s: run time = %.3lf sec (%.1lf min), Q(RF time) = %.5lf C, Q(103.9 kHz time) = %.5lf",
                           var.c_str(),deltaTime,deltaTime/60.,chargeSum,chargeSum103) << std::endl;
      // fill output
      Q.push_back(chargeSum);  
      dQ.push_back(0);
      Time.push_back(deltaTime); 
      dTime.push_back(0); 
      // clear for next run 
      runData.clear();
      chargeSum=0;
      chargeSum103=0;
   }

   return 0;
}
//______________________________________________________________________________
int GetCharge(std::vector<std::string> var,BCMManager *mgr,
              std::vector<double> &run,std::vector<std::vector<double>> &mean,std::vector<std::vector<double>> &stdev){
   // get the charge associated with each run
  
   std::vector<int> rr;
   mgr->GetRunList(rr); 

   // // remove duplicates if necessary 
   // auto last = std::unique(rr.begin(), rr.end());
   // rr.erase(last,rr.end());

   const int NV = var.size();
   double chargeSum[NV],chargeSum103[NV];  
   double MICROAMPS = 1E-6;  

   int M=0;
   double startTime=0,endTime=0,deltaTime=0;
   double timeStep=0,timeStep103=0;
   std::vector<scalerData_t> runData; 
   const int NR = rr.size();
   for(int i=0;i<NR;i++){
      // get data for the run
      mgr->GetVector_scaler("sbs",rr[i],runData);
      // get start and end time of the run
      M = runData.size();
      startTime = runData[0].time; 
      endTime   = runData[M-1].time;
      deltaTime = endTime - startTime;
      // accumulate charge over the whole run
      // for each variable  
      for(int j=0;j<M;j++){
         if(j==0){
	    timeStep    = runData[1].time - runData[0].time;  
	    timeStep103 = runData[1].time - runData[0].time;  
	 }else{
	    timeStep    = runData[j].time - runData[j-1].time;
	    timeStep103 = runData[j].time - runData[j-1].time;
	 }
	 if(j<10) std::cout << Form("event %d, timeStep = %.3lf sec, timeStep(103kHz) = %.3lf",j,timeStep,timeStep103) << std::endl;
         for(int k=0;k<NV;k++){
	    chargeSum[k]    += timeStep*runData[j].getValue(var[k])*MICROAMPS; // convert to Amps 
	    chargeSum103[k] += timeStep103*runData[j].getValue(var[k])*MICROAMPS; // convert to Amps 
         }
      }
      // print to screen 
      std::cout << Form("run %d:",rr[i]) << std::endl;
      for(int k=0;k<NV;k++){
	 std::cout << Form("   %15s: run time = %.3lf sec, Q(RF time) = %.5lf C, Q(103.9 kHz time) = %.5lf",
                           var[k].c_str(),deltaTime,chargeSum[k],chargeSum103[k]) << std::endl;
      }
      // clear for next run 
      runData.clear();
      for(int k=0;k<NV;k++){
	 chargeSum[k] = 0;
	 chargeSum103[k] = 0;
      } 
   } 
   
   return 0;
}

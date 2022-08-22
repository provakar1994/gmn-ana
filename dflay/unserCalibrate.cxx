// Script for use with Unser calibration data
// Compute the (Beam-on) - (beam-off) values 
// then fit corrected rates vs current to a line 

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"

#include "./include/calibCoeff.h"
#include "./include/codaRun.h"
#include "./src/ABA.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/JSONManager.cxx"

double myFitFunc(double *x,double *p); 

TGraphErrors *GetResiduals(TF1 *fit,std::vector<double> x,std::vector<double> y,std::vector<double> ey); 

int unserCalibrate(){

   int rc=0;
 
   JSONManager *jmgr = new JSONManager("./input/json/unser-calib.json"); 
   std::string datapath   = jmgr->GetValueFromKey_str("datapath"); 
   std::string outpath_cc = jmgr->GetValueFromKey_str("outpath_cc"); 
   double fitMin = jmgr->GetValueFromKey<double>("fitMin"); 
   double fitMax = jmgr->GetValueFromKey<double>("fitMax");
   delete jmgr;  

   std::vector<producedVariable_t> data;
   rc = bcm_util::LoadProducedVariables(datapath.c_str(),data);
   if(rc!=0) return 1;

   std::vector<double> w,aba,abaErr;
   std::vector<double> timeOn,on,onErr; 
   std::vector<double> timeOff,off,offErr; 
   std::vector<double> I,uc,ucErr;
  
   ABA *myABA = new ABA();
   myABA->UseTimeWeight(); 
   myABA->SetVerbosity(1);  

   double arg=0,argErr=0,mean=0,err=0,stdev=0;
 
   // compute (beam-on) - (beam-off) differences
   int M=0;
   int grp_prev = data[0].group; // effective beam current  
   const int N = data.size(); 
   for(int i=0;i<N;i++){
      std::cout << Form("==== GROUP %d ====",data[i].group) << std::endl;
      // check the group 
      if(data[i].group==grp_prev){
	 // group match! store data based on beam state 
	 if(data[i].beam_state.compare("on")==0){
	    timeOn.push_back( data[i].time ); 
	    on.push_back( data[i].mean ); 
	    onErr.push_back( data[i].stdev ); 
         }else if(data[i].beam_state.compare("off")==0){
	    timeOff.push_back( data[i].time ); 
	    off.push_back( data[i].mean ); 
	    offErr.push_back( data[i].stdev ); 
         } 
      }else{
	 // compute ABA stats
         // check 
         M = timeOff.size();
         for(int j=0;j<M;j++){
	    std::cout << Form("off: %.3lf, %.3lf ± %.3lf; on: %.3lf, %.3lf ± %.3lf",
                              timeOff[j],off[j],offErr[j],timeOn[j],on[j],onErr[j]) << std::endl;
         }
	 if(M==1){
	    std::cout << "**** ONLY ONE CYCLE! GROUP " << grp_prev << std::endl;
	    // account for one cycle
	    mean = on[0] - off[0]; 
	    err  = TMath::Sqrt( onErr[0]*onErr[0] + offErr[0]*offErr[0] );  
	 }else{  
	    // compute ABA stats
	    myABA->GetDifference(timeOff,off,offErr,timeOn,on,onErr,aba,abaErr);
	    M = aba.size();
	    for(int j=0;j<M;j++){
	       argErr = abaErr[j]*abaErr[j]; 
	       if(abaErr[j]!=0){
		  w.push_back(1./argErr);
	       }else{
		  w.push_back(1); 
	       }
	    }
	    // compute the weighted mean
	    math_df::GetWeightedMean<double>(aba,w,mean,err);
	    stdev = math_df::GetStandardDeviation<double>(aba); 
	    // we compute (A-B), but we actually want B-A
	    mean *= -1; 
	 } 
	 // store results
	 I.push_back( (double)grp_prev ); 
	 uc.push_back(mean);  
	 ucErr.push_back(err); 
	 std::cout << Form("I = %.1lf uA, u_corr = (%.3lf ± %.3lf) Hz",(double)grp_prev,mean,err) << std::endl;
	 std::cout << Form("-----------------------------------------") << std::endl;

	 // set up for next data set 
	 w.clear();
	 timeOn.clear(); 
	 on.clear(); 
	 onErr.clear(); 
	 timeOff.clear(); 
	 off.clear(); 
	 offErr.clear(); 
	 aba.clear(); 
	 abaErr.clear();
	 // store this one since it's needed for the next set!  
	 if(data[i].beam_state.compare("on")==0){
	    timeOn.push_back( data[i].time ); 
	    on.push_back( data[i].mean ); 
	    onErr.push_back( data[i].stdev ); 
         }else if(data[i].beam_state.compare("off")==0){
	    timeOff.push_back( data[i].time ); 
	    off.push_back( data[i].mean ); 
	    offErr.push_back( data[i].stdev ); 
         } 
      }
      grp_prev = data[i].group;
   }
  
   M = timeOff.size();
   for(int j=0;j<M;j++){
      std::cout << Form("off: %.3lf, %.3lf ± %.3lf; on: %.3lf, %.3lf ± %.3lf",
	    timeOff[j],off[j],offErr[j],timeOn[j],on[j],onErr[j]) << std::endl;
   } 
   // new group; compute ABA stats
   w.clear();
   aba.clear(); 
   abaErr.clear();
   if(M==1){
      std::cout << "**** ONLY ONE CYCLE! GROUP " << grp_prev << std::endl;
      // account for one cycle
      mean = on[0] - off[0]; 
      err  = TMath::Sqrt( onErr[0]*onErr[0] + offErr[0]*offErr[0] );  
   }else{  
      // compute ABA stats
      myABA->GetDifference(timeOff,off,offErr,timeOn,on,onErr,aba,abaErr);
      M = aba.size();
      for(int j=0;j<M;j++){
	 argErr = abaErr[j]*abaErr[j]; 
	 if(abaErr[j]!=0){
	    w.push_back(1./argErr);
	 }else{
	    w.push_back(1); 
	 }
      }
      // compute the weighted mean
      math_df::GetWeightedMean<double>(aba,w,mean,err);
      stdev = math_df::GetStandardDeviation<double>(aba); 
      // we compute (A-B), but we actually want B-A
      mean *= -1; 
   }
   // store results
   I.push_back( (double)grp_prev ); 
   uc.push_back(mean); 
   ucErr.push_back(err); 
   std::cout << Form("I = %.1lf uA, u_corr = (%.3lf ± %.3lf) Hz",(double)grp_prev,mean,err) << std::endl;
   std::cout << Form("-----------------------------------------") << std::endl; 
 
   // now make a plot and fit 
   TGraphErrors *g = graph_df::GetTGraphErrors(I,uc,ucErr); 
   graph_df::SetParameters(g,20,kBlack); 

   // set up fit function 
   const int npar = 2;
   double min = fitMin;
   double max = fitMax; // in uA   
   TF1 *myFit = new TF1("myFit",myFitFunc,min,max,npar);
   for(int i=0;i<npar;i++) myFit->SetParameter(i,0);

   TString Title      = Form("Unser Calibration");
   TString xAxisTitle = Form("Beam Current [#muA]");
   TString yAxisTitle = Form("Corrected Unser Rate [Hz]");

   TCanvas *c1 = new TCanvas("c1","Unser Calibration",1200,800);
   c1->Divide(1,2);
   
   c1->cd(1);
   gStyle->SetOptFit(111);
   g->Draw("ap");
   graph_df::SetLabels(g,Title,xAxisTitle,yAxisTitle); 
   graph_df::SetLabelSizes(g,0.05,0.06); 
   g->Draw("ap");
   g->Fit("myFit","QR"); 
   c1->Update();

   TGraphErrors *gr = GetResiduals(myFit,I,uc,ucErr); 
   graph_df::SetParameters(gr,20,kBlack); 
   
   c1->cd(2);
   gr->Draw("alp");
   graph_df::SetLabels(gr,"Fit Residuals",xAxisTitle,"data - fit [Hz]"); 
   graph_df::SetLabelSizes(gr,0.05,0.06); 
   gr->Draw("alp");
   c1->Update();

   // get fit results
   // double arg=0,argErr=0;
   // std::vector<double> par,parErr; 
   calibCoeff_t apt; 
   for(int i=0;i<npar;i++){
      arg    = myFit->GetParameter(i); 
      argErr = myFit->GetParError(i);
      apt.dev = "unser"; 
      if(i==0){
	 apt.offset    = arg; 
	 apt.offsetErr = argErr;
      }else if(i==1){
	 apt.slope    = arg; 
	 apt.slopeErr = argErr;
      }
      std::cout << Form("p[%d] = %.3lf ± %.3lf",i,arg,argErr) << std::endl; 
   }

   std::vector<calibCoeff_t> cc; 
   cc.push_back(apt); 

   bcm_util::WriteToFile_cc(outpath_cc.c_str(),cc);  

   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetResiduals(TF1 *fit,std::vector<double> x,std::vector<double> y,std::vector<double> ey){
   double arg=0,arg_fit=0;
   std::vector<double> r;
   const int N = x.size();
   for(int i=0;i<N;i++){
      arg_fit = fit->Eval(x[i]);
      arg     = y[i] - arg_fit; 
      r.push_back(arg);  
   }

   TGraphErrors *g = graph_df::GetTGraphErrors(x,r,ey);
   return g;
}
//______________________________________________________________________________
double myFitFunc(double *x,double *p){
   // linear fit
   double f = p[0] + p[1]*x[0];
   return f;
}

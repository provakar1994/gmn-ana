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

#include "./include/codaRun.h"
#include "./src/ABA.cxx"
#include "./src/BCMManager.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/Graph.cxx"

TGraphErrors *GetResiduals(TF1 *fit,std::vector<double> x,std::vector<double> y,std::vector<double> ey); 
double myFitFunc(double *x,double *p); 

int unserCalibrate_simple(){

   int rc=0;
   std::vector<std::string> label,path;
   rc = bcm_util::LoadConfigPaths("./input/unser-calib_simple.csv",label,path);

   std::string datapath,outpath_cc;
   const int NF = label.size();
   for(int i=0;i<NF;i++){
      if(label[i].compare("datapath")==0)   datapath   = path[i];
      if(label[i].compare("outpath_cc")==0) outpath_cc = path[i];
   }

   std::vector<producedVariable_t> data;
   rc = bcm_util::LoadProducedVariables(datapath.c_str(),data);
   if(rc!=0) return 1;

   std::vector<double> w,diff;
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
         // check 
         M = timeOff.size();
         for(int j=0;j<M;j++){
	    std::cout << Form("off: %.3lf, %.3lf ± %.3lf; on: %.3lf, %.3lf ± %.3lf",
                              timeOff[j],off[j],offErr[j],timeOn[j],on[j],onErr[j]) << std::endl;
         } 
	 // new group; compute stats
	 // simple method: just take on-off, no weighting
	 for(int j=0;j<M;j++){
	    arg = on[j] - off[j];
	    argErr = TMath::Sqrt( onErr[j]*onErr[j] + offErr[j]*offErr[j] ); 
	    diff.push_back(arg);
	    w.push_back(1./(argErr*argErr));
         } 
	 // compute the weighted mean
         math_df::GetWeightedMean<double>(diff,w,mean,err);
	 stdev = math_df::GetStandardDeviation<double>(diff); 
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
	 diff.clear(); 
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
   // new group; compute stats
   // simple method: just take on-off, no weighting
   diff.clear();
   for(int j=0;j<M;j++){
      arg = on[j] - off[j];
      argErr = TMath::Sqrt( onErr[j]*onErr[j] + offErr[j]*offErr[j] ); 
      diff.push_back(arg);
      w.push_back(1./(argErr*argErr));
   } 
   // compute the weighted mean
   math_df::GetWeightedMean<double>(diff,w,mean,err);
   stdev = math_df::GetStandardDeviation<double>(diff); 
   // compute the weighted mean
   math_df::GetWeightedMean<double>(diff,w,mean,err);
   stdev = math_df::GetStandardDeviation<double>(diff);  
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
   double min = 0;
   double max = 150; // in uA   
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
   g->Fit("myFit","Q");
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
   std::vector<double> par,parErr;
   for(int i=0;i<npar;i++){
      par.push_back( myFit->GetParameter(i) );
      parErr.push_back( myFit->GetParError(i) );
      std::cout << Form("p[%d] = %.3lf ± %.3lf",i,par[i],parErr[i]) << std::endl;
   }

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

// Compute BCM calibration coefficients
// - Uses the output of bcmCalibProcessCuts.cxx
// - Computes (beam-on) - (beam-off) BCM rates (all variables) 
//   and plots as a function of the Unser current.
//   Linear fit produces the calibration coefficients  

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
#include "./src/Utilities.cxx"
#include "./src/bcmUtilities.cxx"

std::string LOG_PATH="";

double myFitFunc(double *x,double *p); 
TGraphErrors *GetTGraphErrors(std::vector<producedVariable_t> unser,std::vector<producedVariable_t> bcm); 

int bcmCalibrate(const char *confPath){

   int rc=0;

   // load configuration file 
   JSONManager *jmgr = new JSONManager(confPath);
   std::string unsccPath = jmgr->GetValueFromKey_str("unser-cc-path"); 
   std::string tag       = jmgr->GetValueFromKey_str("tag"); 

   std::vector<double> fitMin,fitMax; 
   jmgr->GetVectorFromKey<double>("fit-min",fitMin); 
   jmgr->GetVectorFromKey<double>("fit-max",fitMax); 
   delete jmgr;

   // set up output directory paths
   char data_dir[200],plot_dir[200],log_path[200],plt_path[200],outpath[200];
   sprintf(data_dir,"./output/%s"                  ,tag.c_str());
   sprintf(log_path,"./output/%s/log/calibrate.txt",tag.c_str()); 
   sprintf(plot_dir,"./output/%s/plots"            ,tag.c_str());
   LOG_PATH = log_path;
   
   std::string PLOT_PATH=""; 
   
   util_df::LogMessage(log_path,"========== COMPUTING BCM CALIBRATION COEFFICIENTS ==========",'a');

   // load Unser data 
   char unserPath[200]; 
   sprintf(unserPath,"%s/unser.csv",data_dir); 
   std::vector<producedVariable_t> unser;
   rc = bcm_util::LoadProducedVariables(unserPath,unser);
   if(rc!=0) return 1;

   // compute pedestal rates
   double mean_ped=0,stdev_ped=0;
   std::vector<producedVariable_t> unser_off;
   bcm_util::CalculateStatsForBeamState("off",unser,unser_off,mean_ped,stdev_ped,LOG_PATH);
   sprintf(outpath,"%s/unser_ped.csv",data_dir);
   bcm_util::WriteToFile(outpath,unser_off);  

   // compute pedestal-subtracted Unser rates and convert to current
   std::vector<producedVariable_t> unser_ps;
   sprintf(plt_path,"%s/unser.pdf",plot_dir); 
   PLOT_PATH = plt_path; 
   bcm_util::CalculatePedestalSubtraction(unser,unser_ps,LOG_PATH,PLOT_PATH); 
   // load Unser calibration coefficients 
   std::vector<calibCoeff_t> unsCC;
   bcm_util::LoadFittedOffsetGainData(unsccPath.c_str(),unsCC); 
   // convert to current 
   std::vector<producedVariable_t> unser_cur;
   bcm_util::ConvertToCurrent(unsCC[0],unser_ps,unser_cur); // only ONE set of unser calibration coefficients!  

   // load each BCM variable; compute pedestal-subtracted 
   // values. create plots of ped-subtracted rates vs Unser current 
   const int N = 6; 
   std::string bcmVar[N] = {"u1","unew","d1","d3","d10","dnew"};

   TGraphErrors **g = new TGraphErrors*[N]; 
   TF1 **myFit      = new TF1*[N]; 
 
   // set up fit parameters
   const int npar=2;
   double offset=0,offsetErr=0,slope=0,slopeErr=0; 

   // calibration coefficient output 
   calibCoeff_t ccPt; 
   std::vector<calibCoeff_t> cc; 

   int NCOL=3,NROW=2;
   TCanvas *c1 = new TCanvas("c1","BCM Calibration",1200,800);
   c1->Divide(NCOL,NROW);  

   TString Title,xAxisTitle,yAxisTitle,fitName;

   char inpath[200],msg[200]; 
   
   std::vector<producedVariable_t> bcm,bcm_ps,bcm_off; 
   for(int i=0;i<N;i++){
      // load data
      sprintf(inpath,"%s/%s.csv",data_dir,bcmVar[i].c_str()); 
      rc = bcm_util::LoadProducedVariables(inpath,bcm);
      // compute pedestal rates
      bcm_util::CalculateStatsForBeamState("off",bcm,bcm_off,mean_ped,stdev_ped,LOG_PATH);
      sprintf(outpath,"%s/%s_ped.csv",data_dir,bcmVar[i].c_str());
      bcm_util::WriteToFile(outpath,bcm_off);  
      // subtract pedestal
      sprintf(plt_path,"%s/%s.pdf",plot_dir,bcmVar[i].c_str()); 
      PLOT_PATH = plt_path;
      bcm_util::CalculatePedestalSubtraction(bcm,bcm_ps,LOG_PATH,PLOT_PATH);
      // create TGraphError plot  
      g[i] = GetTGraphErrors(unser_cur,bcm_ps); 
      graph_df::SetParameters(g[i],20,kBlack);
      // set up the fit function 
      fitName  = Form("fit_%s",bcmVar[i].c_str());
      myFit[i] = new TF1(fitName,myFitFunc,fitMin[i],fitMax[i],npar);
      for(int j=0;j<npar;j++){
	 myFit[i]->SetParameter(j,0);
         myFit[i]->SetParLimits(j,-10E+3,10E+3);
      }
      // set titles 
      Title      = Form("%s",bcmVar[i].c_str());
      xAxisTitle = Form("Unser Current [#muA]");
      yAxisTitle = Form("%s rate [Hz]",bcmVar[i].c_str());
      // plot data 
      c1->cd(i+1);  
      gStyle->SetOptFit(111); 
      g[i]->Draw("ap"); 
      graph_df::SetLabels(g[i],Title,xAxisTitle,yAxisTitle);
      g[i]->Draw("ap");
      // fit the data
      g[i]->Fit(fitName,"QR");
      c1->Update(); 
      // extract fit results
      offset    = myFit[i]->GetParameter(0); 
      offsetErr = myFit[i]->GetParError(0); 
      slope     = myFit[i]->GetParameter(1); 
      slopeErr  = myFit[i]->GetParError(1); 
      // store results 
      ccPt.dev         = bcmVar[i];
      ccPt.pedestal    = mean_ped;  // averaged over all cycles and groups  
      ccPt.pedestalErr = stdev_ped;
      ccPt.offset      = offset; 
      ccPt.offsetErr   = offsetErr;
      ccPt.slope       = slope;
      ccPt.slopeErr    = slopeErr;
      cc.push_back(ccPt);
      // print to screen 
      sprintf(msg,"bcm = %s, pedestal = %.3E ± %.3E, offset = %.3E ± %.3E, slope = %.3E ± %.3E",
                        ccPt.dev.c_str(),ccPt.pedestal,ccPt.pedestalErr,ccPt.offset,ccPt.offsetErr,ccPt.slope,ccPt.slopeErr);
      util_df::LogMessage(log_path,msg,'a'); 
      // set up for next BCM
      bcm.clear();
      bcm_ps.clear(); 
      bcm_off.clear(); 
   }

   util_df::LogMessage(log_path,"----------------",'a');

   // save the canvas
   TString plotPath = Form("%s/%s.pdf",plot_dir,tag.c_str()); 
   c1->cd();
   c1->Print(plotPath);  

   // print results to file
   sprintf(outpath,"%s/result.csv",data_dir); 
   bcm_util::WriteToFile_cc(outpath,cc);  

   return 0;
}
//______________________________________________________________________________
TGraphErrors *GetTGraphErrors(std::vector<producedVariable_t> unser,std::vector<producedVariable_t> bcm){
   // produce a plot of the BCM rate (y axis) against the unser current (x axis)  

   const int NX = unser.size();
   const int NY = bcm.size();
   if(NX!=NY){
      util_df::LogMessage(LOG_PATH.c_str(),"[GetTGraphErrors]: ERROR! Unser NPTS != BCM NPTS!",'a');
      exit(1); 
   }

   std::vector<double> x,ex,y,ey; 
   for(int i=0;i<NX;i++){
      x.push_back(unser[i].mean); 
      ex.push_back(unser[i].stdev); 
      y.push_back(bcm[i].mean); 
      ey.push_back(bcm[i].stdev); 
   }

   TGraphErrors *g = graph_df::GetTGraphErrors(x,ex,y,ey); 
   return g;
}
//______________________________________________________________________________
double myFitFunc(double *x,double *p){
   // linear fit
   double f = p[0] + p[1]*x[0];
   return f;
}

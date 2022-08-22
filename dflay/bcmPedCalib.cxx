// Compute the mean value of the BCMs and Unser 
// for beam-off (pedestal) run sets  
// output is pedestal files for each BCM and Unser  

#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"

#include "./include/codaRun.h"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

int GetStats(std::vector<double> x,std::vector<double> dx,double &mean,double &stdev); 
int WriteToFile(const char *outpath,const char *dev,int runMin,int runMax,double x,double dx); 
int WriteToFile_unser(const char *outpath,std::vector<int> runMin,std::vector<int> runMax,
                      std::vector<double> offset,std::vector<double> offsetErr); 

int bcmPedCalib(const char *runPath){

   char oprefix[200]; 
   sprintf(oprefix,"./output/pedestal");

   // settings 
   bool logScale = false;

   gStyle->SetOptStat(0);

   int rc=0;

   TString prefix;
   std::vector<codaRun_t> runList;  
   rc = bcm_util::LoadRuns(runPath,prefix,runList);
   if(rc!=0) return 1; 

   BCMManager *mgr = new BCMManager("NONE",false,"NONE");

   TString filePath;  
   const int NR = runList.size();  
   for(int i=0;i<NR;i++){ 
      filePath = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d.root",
                      prefix.Data(),runList[i].runNumber,runList[i].stream,runList[i].segmentBegin,runList[i].segmentEnd);
      mgr->LoadFile(filePath,runList[i].runNumber);
   }

   // do stats by run number 
   std::vector<scalerData_t> rawData,data;
   mgr->GetVector_scaler("sbs",rawData);   

   // load cuts 
   std::vector<cut_t> cutList; 
   cut_util::LoadCuts("./input/cut-list_bcm-ped.csv",cutList);
   // apply cuts  
   cut_util::ApplyCuts(cutList,rawData,data);

   // sort the data by run number 
   std::sort(data.begin(),data.end(),compareScalerData_byRun);

   const int N = 7; 
   TString var[N] = {"unser.rate","u1.rate","unew.rate","d1.rate","d3.rate","d10.rate","dnew.rate"};

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
      std::cout << Form("%05d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
	                (int)run[i],mean_uns[i],stdev_uns[i],mean_u1[i],stdev_u1[i],
                        mean_unew[i],stdev_unew[i],mean_d1[i],stdev_d1[i],mean_d3[i],stdev_d3[i],
                        mean_d10[i],stdev_d10[i],mean_dnew[i],stdev_dnew[i]) << std::endl;
   }

   // BCM analysis: get mean and std dev across all runs
   // (these are stable) 
   double MEAN_U1=0,STDEV_U1=0,MEAN_UNEW=0,STDEV_UNEW=0; 
   double MEAN_D1=0,STDEV_D1=0,MEAN_DNEW=0,STDEV_DNEW=0; 
   double MEAN_D3=0,STDEV_D3=0,MEAN_D10=0,STDEV_D10=0; 
   GetStats(mean_u1  ,stdev_u1  ,MEAN_U1  ,STDEV_U1);
   GetStats(mean_unew,stdev_unew,MEAN_UNEW,STDEV_UNEW);
   GetStats(mean_d1  ,stdev_d1  ,MEAN_D1  ,STDEV_D1);
   GetStats(mean_d3  ,stdev_d3  ,MEAN_D3  ,STDEV_D3);
   GetStats(mean_d10 ,stdev_d10 ,MEAN_D10 ,STDEV_D10);
   GetStats(mean_dnew,stdev_dnew,MEAN_DNEW,STDEV_DNEW);

   int runMin=0,runMax=99999;

   char outpath_u1[200],outpath_unew[200];
   char outpath_d1[200],outpath_d3[200],outpath_d10[200],outpath_dnew[200];
   sprintf(outpath_u1  ,"%s/pedestal_u1.csv"  ,oprefix); 
   sprintf(outpath_unew,"%s/pedestal_unew.csv",oprefix); 
   sprintf(outpath_d1  ,"%s/pedestal_d1.csv"  ,oprefix); 
   sprintf(outpath_d3  ,"%s/pedestal_d3.csv"  ,oprefix); 
   sprintf(outpath_d10 ,"%s/pedestal_d10.csv" ,oprefix); 
   sprintf(outpath_dnew,"%s/pedestal_dnew.csv",oprefix); 

   WriteToFile(outpath_u1  ,"u1"  ,runMin,runMax,MEAN_U1  ,STDEV_U1  ); 
   WriteToFile(outpath_unew,"unew",runMin,runMax,MEAN_UNEW,STDEV_UNEW); 
   WriteToFile(outpath_d1  ,"d1"  ,runMin,runMax,MEAN_D1  ,STDEV_D1  ); 
   WriteToFile(outpath_d3  ,"d3"  ,runMin,runMax,MEAN_D3  ,STDEV_D3  ); 
   WriteToFile(outpath_d10 ,"d10" ,runMin,runMax,MEAN_D10 ,STDEV_D10 ); 
   WriteToFile(outpath_dnew,"dnew",runMin,runMax,MEAN_DNEW,STDEV_DNEW); 

   // Unser analysis: read in the Unser run groups
   CSVManager *csv = new CSVManager(); 
   rc = csv->ReadFile("./input/unser-ped-calib-grps.csv",true);
   if(rc!=0) return 1; 
   
   std::vector<int> RUN_MIN,RUN_MAX; 
   csv->GetColumn_byName<int>("runMin",RUN_MIN);  
   csv->GetColumn_byName<int>("runMax",RUN_MAX);  

   // loop over all valid runs, identify the run group, store data, get stats
   double mean=0,stdev=0;  
   std::vector<double> rr,x,dx,MU,SIG;  
   int M = RUN_MIN.size();
   int NS=0;
   for(int i=0;i<M;i++){
      for(int j=0;j<NNR;j++){
	 if(run[j]>=RUN_MIN[i]&&run[j]<=RUN_MAX[i]){
	    if(stdev_uns[j]!=0){
	       x.push_back( mean_uns[j] );
	       dx.push_back( stdev_uns[j] ); 
	    } 
	 }
      }
      // gathered all runs for this group, get stats
      NS = x.size();
      GetStats(x,dx,mean,stdev);
      MU.push_back(mean); 
      if(NS==1){
	 // only one run, use stdev of that run
	 std::cout << "WARNING: only 1 run for group " << i+1 << std::endl;  
	 SIG.push_back(dx[0]); 
      }else{
	 SIG.push_back(stdev);
      } 
      std::cout << Form("Group %d: %d--%d, %.3lf, %.3lf",i+1,RUN_MIN[i],RUN_MAX[i],MU[i],SIG[i]) << std::endl;  
      // set up for next group
      rr.clear();
      x.clear();
      dx.clear();
   }

   char outpath_uns[200]; 
   sprintf(outpath_uns,"%s/pedestal_unser.csv",oprefix); 
   WriteToFile_unser(outpath_uns,RUN_MIN,RUN_MAX,MU,SIG);

   return 0;
}
//______________________________________________________________________________
int GetStats(std::vector<double> x,std::vector<double> dx,double &mean,double &stdev){
   // get weighted avg and standard deviation 

   double arg=0;  
   const int N = x.size(); 
   std::vector<double> w; 
   for(int i=0;i<N;i++){
      arg = 1; 
      if(dx[i]!=0) arg = 1./(dx[i]*dx[i]);
      w.push_back(arg);  
   } 

   int rc=0; 
   double err=0; 
   rc    = math_df::GetWeightedMean<double>(x,w,mean,err);
   stdev = math_df::GetStandardDeviation<double>(x);  

   return 0;
}
//______________________________________________________________________________
int WriteToFile(const char *outpath,const char *dev,int runMin,int runMax,double x,double dx){
   // for BCMs (not Unser) 
   std::string header = "dev,runMin,runMax,pedestal,pedestalErr";

   std::vector<std::string> DEV;
   std::vector<int> RUN_MIN,RUN_MAX; 
   std::vector<double> pedestal,pedestalErr,offset,offsetErr,gain,gainErr;
   pedestal.push_back(x); 
   pedestalErr.push_back(dx);
   RUN_MIN.push_back(runMin);   
   RUN_MAX.push_back(runMax);  
   DEV.push_back(dev);  
 
   int NROW = RUN_MIN.size(); 
   int NCOL = 5; 

   CSVManager *csv = new CSVManager();
   csv->InitTable(NROW,NCOL);
   csv->SetColumn_str(0,DEV);
   csv->SetColumn<int>(1,RUN_MIN); 
   csv->SetColumn<int>(2,RUN_MAX); 
   csv->SetColumn<double>(3,pedestal);
   csv->SetColumn<double>(4,pedestalErr);
   csv->SetHeader(header);
   csv->WriteFile(outpath);

   delete csv;
   return 0;
}
//______________________________________________________________________________
int WriteToFile_unser(const char *outpath,std::vector<int> runMin,std::vector<int> runMax,
                      std::vector<double> pedestal,std::vector<double> pedestalErr){
   // for Unser 
   std::string header = "dev,runMin,runMax,pedestal,pedestalErr";

   std::vector<std::string> dev; 
   
   const int NROW = runMin.size();
   int NCOL = 5; 

   for(int i=0;i<NROW;i++){
      dev.push_back("unser");  
   }

   CSVManager *csv = new CSVManager();
   csv->InitTable(NROW,NCOL);
   csv->SetColumn_str(0,dev);
   csv->SetColumn<int>(1,runMin); 
   csv->SetColumn<int>(2,runMax); 
   csv->SetColumn<double>(3,pedestal);
   csv->SetColumn<double>(4,pedestalErr);
   csv->SetHeader(header);
   csv->WriteFile(outpath);

   delete csv;
   return 0;
}

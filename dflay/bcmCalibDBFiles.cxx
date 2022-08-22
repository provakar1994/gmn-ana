// read in pedestal and fitted offset and slope (gain) files 
// and create calibration coefficient files
// for use in the analyzer 

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

int GetStats_pv(std::vector<producedVariable_t> data,double &mean,double &err,double &stdev); 

int bcmCalibDBFiles(const char *confPath){
   
   int rc=0;

   JSONManager *jmgr = new JSONManager(confPath);
   std::string tag        = jmgr->GetValueFromKey_str("tag");
   std::string unser_path = jmgr->GetValueFromKey_str("unser-cc-path");   
   delete jmgr;  
 
   char inpath[200],prefix[200],log_path[200],db_path_lhrs[200],db_path_sbs[200],msg[200]; 
   sprintf(prefix,"./output/%s",tag.c_str());
   sprintf(db_path_lhrs,"%s/db_LeftBCM.dat",prefix);
   sprintf(db_path_sbs ,"%s/db_sbsBCM.dat",prefix);
   sprintf(log_path    ,"%s/log/db.txt"   ,prefix);

   // load ALL results 
   std::vector<calibCoeff_t> cc;
   sprintf(inpath,"%s/result.csv",prefix);  
   rc = bcm_util::LoadCalibrationCoefficients(inpath,cc); 

   // also get the Unser
   // precision current calibration 
   rc = bcm_util::LoadCalibrationCoefficients(unser_path.c_str(),cc);
   const int NCC = cc.size();  

   // unser pedestal data from this analysis 
   sprintf(inpath,"%s/unser_ped.csv",prefix);
   std::vector<producedVariable_t> uped; 
   rc = bcm_util::LoadProducedVariables(inpath,uped);
   // averaged pedestal 
   double mean=0,err=0,stdev=0;
   GetStats_pv(uped,mean,err,stdev);
   // set unser pedestal 
   cc[NCC-1].pedestal    = mean; 
   cc[NCC-1].pedestalErr = stdev; 

   char dev_msg[200]; 
   char gain_msg[200],offs_msg[200],satq_msg[200];
   char sato_msg[200],curt_msg[200],curi_msg[200]; 

   double totOffset=0,totOffsetErr=0;
   for(int i=0;i<NCC;i++){
      // DB expects a TOTAL offset -- that is pedestal + fitted offset
      // we construct it here 
      totOffset    = cc[i].pedestal + cc[i].offset;
      totOffsetErr = TMath::Sqrt( cc[i].pedestalErr*cc[i].pedestalErr + cc[i].offsetErr*cc[i].offsetErr ); 
      // sprintf(msg,"%s: pedestal = %.3lf, offset = %.3lf, total = %.3lf",cc[i].dev.c_str(),cc[i].pedestal,cc[i].offset,totOffset); 
      // util_df::LogMessage(log_path,msg,'a'); 
      if(i==0){
	 sprintf(dev_msg ,"BCM_Names                   = %s",cc[i].dev.c_str()); 
	 sprintf(gain_msg,"BCM_Gain                    = %.3lf",cc[i].slope); 
	 sprintf(offs_msg,"BCM_Offset                  = %.3lf",totOffset); 
	 sprintf(satq_msg,"BCM_SatQuadratic            = %.3lf",0.); 
	 sprintf(sato_msg,"BCM_SatOffset               = %.3lf",0.); 
	 sprintf(curt_msg,"BCM_Current_threshold       = %.3lf",0.); 
	 sprintf(curi_msg,"BCM_Current_threshold_index = %d"   ,0 ); 
      }else{
	 sprintf(dev_msg ,"%s %s"   ,dev_msg ,cc[i].dev.c_str()); 
	 sprintf(gain_msg,"%s %.3lf",gain_msg,cc[i].slope); 
	 sprintf(offs_msg,"%s %.3lf",offs_msg,totOffset); 
	 sprintf(satq_msg,"%s %.3lf",satq_msg,0.); 
	 sprintf(sato_msg,"%s %.3lf",sato_msg,0.); 
	 sprintf(curt_msg,"%s %.3lf",curt_msg,0.); 
	 sprintf(curi_msg,"%s %d"   ,curi_msg,0 ); 
      }
   }

   char nbcm_msg[200]; 
   sprintf(nbcm_msg,"NumBCMs                     = %d",NCC); 

   std::string nbcm_msg_str = nbcm_msg;
   std::string dev_msg_str  = dev_msg;  
   std::string gain_msg_str = gain_msg;
   std::string offs_msg_str = offs_msg;
   std::string satq_msg_str = satq_msg;
   std::string sato_msg_str = sato_msg;
   std::string curt_msg_str = curt_msg;     
   std::string curi_msg_str = curi_msg;

   std::vector<std::string> MSG;
   MSG.push_back(nbcm_msg_str);  
   MSG.push_back(dev_msg_str); 
   MSG.push_back(gain_msg_str); 
   MSG.push_back(offs_msg_str); 
   MSG.push_back(satq_msg_str); 
   MSG.push_back(sato_msg_str); 
   MSG.push_back(curt_msg_str); 
   MSG.push_back(curi_msg_str); 

   const int NS = MSG.size();
   for(int i=0;i<NS;i++){
      util_df::LogMessage(db_path_lhrs,MSG[i].c_str(),'a'); 
      util_df::LogMessage(db_path_sbs ,MSG[i].c_str(),'a'); 
   }

   return 0;
}
//______________________________________________________________________________
int GetStats_pv(std::vector<producedVariable_t> data,double &mean,double &err,double &stdev){
   // get the mean, std-err of mean, stdev on a vector of producedVariable data types
   double argErr=0;
   const int N = data.size();
   std::vector<double> v,w; 
   for(int i=0;i<N;i++){
      argErr = data[i].stdev; 
      if( argErr!=0 ){
	 w.push_back( 1./(argErr*argErr) ); 
      }else{
	 w.push_back(1); 
      }
      v.push_back(data[i].mean);
   } 
   // calculate stats
   math_df::GetWeightedMean<double>(v,w,mean,err); 
   stdev = math_df::GetStandardDeviation<double>(v); 
   return 0;
}

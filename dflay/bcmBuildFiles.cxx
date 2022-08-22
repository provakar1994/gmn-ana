// read in pedestal and fitted offset and slope (gain) files 
// and create calibration coefficient files 

#include <cstdlib>
#include <iostream>

#include "./src/Graph.cxx"
#include "./src/JSONManager.cxx"
#include "./src/bcmUtilities.cxx"

int bcmBuildFiles(const char *confPath){
   
   int rc=0;

   JSONManager *jmgr = new JSONManager(confPath);
   std::string inpath_unser = jmgr->GetValueFromKey_str("unser-path"); 
   std::string inpath_bcm   = jmgr->GetValueFromKey_str("bcm-path");
   delete jmgr;  

   const int N = 7; 
   std::string bcm[N] = {"unser","u1","unew","d1","d3","d10","dnew"};
 
   char prefix_ped[200],prefix_cc[200]; 
   sprintf(prefix_ped,"./output/pedestal");
   sprintf(prefix_cc ,"./output/calib-coeff");

   std::vector<calibCoeff_t> ped,og,cc; 

   char inpath[200],outpath[200]; 

   int M=0;
   for(int i=0;i<N;i++){
      // get pedestal data for a given BCM 
      sprintf(inpath ,"%s/pedestal_%s.csv"   ,prefix_ped,bcm[i].c_str()); 
      rc = bcm_util::LoadPedestalData(inpath,ped);
      // build the *full* calibration coefficient struct
      // first we copy all pedestal info to the new cc vector 
      M = ped.size();  // determine number of pedestal periods 
      for(int j=0;j<M;j++) cc.push_back(ped[j]);
      // now load the fitted offsets and gains
      if(bcm[i].compare("unser")==0){
	 sprintf(inpath,"%s",inpath_unser.c_str());
      }else{
	 sprintf(inpath,"%s/%s-calib-results.csv",inpath_bcm.c_str(),bcm[i].c_str());
      }
      // apply calibration coeffients  
      rc = bcm_util::LoadFittedOffsetGainData(inpath,og); 
      for(int j=0;j<M;j++){
	 cc[j].offset    = og[0].offset;
	 cc[j].offsetErr = og[0].offsetErr;
	 cc[j].slope     = og[0].slope;
	 cc[j].slopeErr  = og[0].slopeErr;
      }
      // print to file 
      sprintf(outpath,"%s/01-18-22/calib-coeff_%s.csv",prefix_cc,bcm[i].c_str()); 
      bcm_util::WriteToFile_cc(outpath,cc);
      // set up for next BCM 
      cc.clear();
      og.clear();
      ped.clear(); 
   }  

   return 0;
}

#include "../include/BCMManager.h"
//______________________________________________________________________________
BCMManager::BCMManager(const char *filePath,const char *ccDirPath,bool isDebug){
   fIsDebug            = isDebug; 
   fCalculateCurrent   = false; 
   fVerbosity          = 0;  
   fEvtCntrLeft        = 0;  
   fEvtCntrSBS         = 0;  
   fEvtCntrEPICS       = 0; 
   fLastTimeLeft       = 0; 
   fLastTimeSBS        = 0;
   fLastRun            = 0; 
   fLastRunEvtCntLeft  = 0; 
   fLastRunEvtCntSBS   = 0; 
   fLastRunEvtCntEPICS = 0;
   fTimeStamp          = "NONE"; 
   fUTCTimeStamp       = 0;  

   std::string path = filePath; 
   if(path.compare("NONE")!=0){
      LoadFile(filePath);
   }

   std::string cc_dirPath = ccDirPath; 
   if(cc_dirPath.compare("NONE")!=0){
      LoadCalibrationCoefficients(ccDirPath);
   } 
 
}
//______________________________________________________________________________
BCMManager::~BCMManager(){
   Clear();
}
//______________________________________________________________________________
void BCMManager::Clear(){
   fLeft.clear();
   fSBS.clear();
   fEPICS.clear();
   fRunList.clear(); 
   fEvtCntrLeft  = 0;
   fEvtCntrSBS   = 0;
   fEvtCntrEPICS = 0;
   fLastTimeLeft = 0; 
   fLastTimeSBS  = 0;
   fTimeStamp    = "NONE";  
}
//______________________________________________________________________________
int BCMManager::LoadFile(const char *filePath,int runNumber){
   // Attach trees for LHRS and SBS
   int fail=0;
   int rc = CheckFile(filePath); 
   if(rc!=0){
      std::cout << "[BCMManager::LoadFile]: Cannot open the file for run " << runNumber << std::endl;
      return rc;
   }

   // file is valid, read in the data from the trees 
   int rc_lhrs  = LoadDataFromTree(filePath,"TSLeft",runNumber);
   int rc_sbs   = LoadDataFromTree(filePath,"TSsbs" ,runNumber);
   int rc_epics = LoadEPICSDataFromTree(filePath,runNumber); 

   // add the run to the run list
   // must at least have the SBS branch to be considered good
   // also make sure we don't have duplicates
   int notFound=0;
   int NR=fRunList.size(); 
   if(rc_sbs==0){
      std::cout << Form("-----------------------------------------") << std::endl;
      std::cout << Form("[BCMManager::LoadFile]: Load successful: ") << std::endl;
      std::cout << Form("Run:  %d"      ,runNumber) << std::endl;
      std::cout << Form("Date: %s (%lu)",fTimeStamp.c_str(),fUTCTimeStamp) << std::endl;
      std::cout << Form("-----------------------------------------") << std::endl;
      // std::cout << Form("[BCMManager::LoadFile]: Loaded run %d, recorded on %s (%lu)",runNumber,fTimeStamp.c_str(),fUTCTimeStamp) << std::endl;
      for(int i=0;i<NR;i++){
	 if(runNumber!=fRunList[i]) notFound += 1;
      }
      if(notFound==NR) fRunList.push_back(runNumber); // all runs don't match the new one, can add it 
   }

   fLastRun = runNumber;

   return rc; 
}
//______________________________________________________________________________
int BCMManager::GetTimeStamp(TFile *f,std::string &timeStamp,unsigned long &utc){

   // get date and time of run
   std::string dayStr="";  
   THaRun* run      = (THaRun*)f->Get("Run_Data");
   TDatime run_date = run->GetDate();
   int dayOfWeek = run_date.GetDayOfWeek();
   int month     = run_date.GetMonth(); 
   int day       = run_date.GetDay();
   int year      = run_date.GetYear();
   int hour      = run_date.GetHour();
   int min       = run_date.GetMinute();
   int sec       = run_date.GetSecond();
   if(dayOfWeek==1) dayStr = "Mon"; 
   if(dayOfWeek==2) dayStr = "Tue"; 
   if(dayOfWeek==3) dayStr = "Wed"; 
   if(dayOfWeek==4) dayStr = "Thu"; 
   if(dayOfWeek==5) dayStr = "Fri"; 
   if(dayOfWeek==6) dayStr = "Sat"; 
   if(dayOfWeek==7) dayStr = "Sun";

   char timeStr[200];
   sprintf(timeStr,"%s %d-%02d-%02d %02d:%02d:%02d",dayStr.c_str(),year,month,day,hour,min,sec);
   timeStamp = timeStr;

   // now compute the UTC time stamp 
   struct tm timeinfo = {0};
   timeinfo.tm_year = year-1900;  // years since 1900 
   timeinfo.tm_mon  = month-1;    // months since january (0--11)  
   timeinfo.tm_mday = day;
   timeinfo.tm_hour = hour;
   timeinfo.tm_min  = min;
   timeinfo.tm_sec  = sec;
   int corr = 0;
   // FIXME: there's got to be a better way to do this... 
   int isDST = timeinfo.tm_isdst; 
   if(isDST) corr = 3600; // daylight savings time, subtract 1 hour  
   time_t timeSinceEpoch = mktime(&timeinfo) - corr;
   utc = (unsigned long)timeSinceEpoch; 

   return 0;
}
//______________________________________________________________________________
int BCMManager::CheckFile(const char *filePath){
   TFile *myFile = new TFile(filePath);

   Long64_t bytesRead = myFile->GetBytesRead();

   int rc=1;
 
   if(bytesRead!=0){
      // file has non-zero size
      rc = 0;
      // get the time stamp of the run 
      GetTimeStamp(myFile,fTimeStamp,fUTCTimeStamp);
   }
 
   return rc;
}
//______________________________________________________________________________
int BCMManager::LoadDataFromTree(const char *filePath,const char *treeName,int runNumber){

   if(fIsDebug) std::cout << "[BCMManager::LoadDataFromTree]: Loading data from tree: " << treeName << std::endl;

   std::string arm; 
   std::string name = treeName;
   if(name.compare("TSLeft")==0) arm = "Left"; 
   if(name.compare("TSsbs")==0)  arm = "sbs"; 

   Long64_t trigEvt=0;
   double trigEvt2; 
 
   double time=0,time_num=0,time_den=0; 
   double time103kHz=0,time103kHz_num=0,time103kHz_den=0; 
   double unser_rate=0,unser_cnt=0,unser_cur=0;
   double u1_rate=0,u1_cnt=0,u1_cur=0;
   double unew_rate=0,unew_cnt=0,unew_cur=0;
   double dnew_rate=0,dnew_cnt=0,dnew_cur=0;
   double d1_rate=0,d1_cnt=0,d1_cur=0;
   double d3_rate=0,d3_cnt=0,d3_cur=0;
   double d10_rate=0,d10_cnt=0,d10_cur=0;

   TChain *ch = nullptr; 
   ch = new TChain(treeName); 
   ch->Add(filePath);
   int NN = ch->GetEntries();

   if(ch==nullptr){
      return 1;
   } 
 
   TTree *aTree = ch->GetTree();
   if(aTree==nullptr) return 1; 

   aTree->SetBranchAddress("evnum"            ,&trigEvt );
   aTree->SetBranchAddress("evNumber"         ,&trigEvt2 );
   aTree->SetBranchAddress(Form("%s.bcm.unser.cnt" ,arm.c_str()),&unser_cnt );
   aTree->SetBranchAddress(Form("%s.bcm.unser.rate",arm.c_str()),&unser_rate);
   aTree->SetBranchAddress(Form("%s.bcm.unser.current",arm.c_str()),&unser_cur);
   aTree->SetBranchAddress(Form("%s.bcm.u1.cnt"    ,arm.c_str()),&u1_cnt    );
   aTree->SetBranchAddress(Form("%s.bcm.u1.rate"   ,arm.c_str()),&u1_rate   );
   aTree->SetBranchAddress(Form("%s.bcm.u1.current",arm.c_str()),&u1_cur    );
   aTree->SetBranchAddress(Form("%s.bcm.unew.cnt"  ,arm.c_str()),&unew_cnt  );
   aTree->SetBranchAddress(Form("%s.bcm.unew.rate" ,arm.c_str()),&unew_rate );
   aTree->SetBranchAddress(Form("%s.bcm.unew.current",arm.c_str()),&unew_cur );
   aTree->SetBranchAddress(Form("%s.bcm.d1.cnt"    ,arm.c_str()),&d1_cnt    );
   aTree->SetBranchAddress(Form("%s.bcm.d1.rate"   ,arm.c_str()),&d1_rate   );
   aTree->SetBranchAddress(Form("%s.bcm.d1.current",arm.c_str()),&d1_cur    );
   aTree->SetBranchAddress(Form("%s.bcm.d3.cnt"    ,arm.c_str()),&d3_cnt    );
   aTree->SetBranchAddress(Form("%s.bcm.d3.rate"   ,arm.c_str()),&d3_rate   );
   aTree->SetBranchAddress(Form("%s.bcm.d3.current",arm.c_str()),&d3_cur    );
   aTree->SetBranchAddress(Form("%s.bcm.d10.cnt"   ,arm.c_str()),&d10_cnt   );
   aTree->SetBranchAddress(Form("%s.bcm.d10.rate"  ,arm.c_str()),&d10_rate  );
   aTree->SetBranchAddress(Form("%s.bcm.d10.current",arm.c_str()),&d10_cur  );
   aTree->SetBranchAddress(Form("%s.bcm.dnew.cnt"  ,arm.c_str()),&dnew_cnt  );
   aTree->SetBranchAddress(Form("%s.bcm.dnew.rate" ,arm.c_str()),&dnew_rate );
   aTree->SetBranchAddress(Form("%s.bcm.dnew.current",arm.c_str()),&dnew_cur );
   if(arm.compare("Left")==0){
      aTree->SetBranchAddress(Form("Left.104kHz_CLK.cnt")       ,&time_num);
      aTree->SetBranchAddress(Form("Left.104kHz_CLK.rate")      ,&time_den);
      aTree->SetBranchAddress(Form("Left.104kHz_CLK.cnt")       ,&time103kHz_num);
      aTree->SetBranchAddress(Form("Left.104kHz_CLK.rate")      ,&time103kHz_den);
   }else if(arm.compare("sbs")==0){
      aTree->SetBranchAddress(Form("sbs.BBCalHi.RF.scaler")     ,&time_num);
      aTree->SetBranchAddress(Form("sbs.BBCalHi.RF.scalerRate") ,&time_den);
      aTree->SetBranchAddress(Form("sbs.104kHz_CLK.cnt")        ,&time103kHz_num);
      aTree->SetBranchAddress(Form("sbs.104kHz_CLK.rate")       ,&time103kHz_den);
   }

   scalerData_t pt; 

   for(int i=0;i<NN;i++){
      aTree->GetEntry(i);
      if(time_den!=0){
	 time = (double)fUTCTimeStamp + time_num/time_den; 
      }else{
	 time = 0;  
      } 
      if(time103kHz_den!=0){
	 time103kHz = (double)fUTCTimeStamp + time103kHz_num/time103kHz_den;
      }else{
	 time103kHz = 0;
      }
      pt.triggerEvent   = trigEvt;
      pt.triggerEvent2  = (signed long long)trigEvt2;
      pt.time_num       = time_num;  
      pt.time_den       = time_den;  
      pt.time103kHz_num = time103kHz_num;  
      pt.time103kHz_den = time103kHz_den;  
      pt.arm            = arm; 
      pt.runNumber      = runNumber; 
      pt.unserRate      = unser_rate;  
      pt.unserCounts    = unser_cnt;  
      pt.u1Rate         = u1_rate;  
      pt.u1Counts       = u1_cnt;  
      pt.unewRate       = unew_rate;  
      pt.unewCounts     = unew_cnt;  
      pt.dnewRate       = dnew_rate;  
      pt.dnewCounts     = dnew_cnt;  
      pt.d1Rate         = d1_rate;  
      pt.d1Counts       = d1_cnt;  
      pt.d3Rate         = d3_rate;  
      pt.d3Counts       = d3_cnt;  
      pt.d10Rate        = d10_rate;  
      pt.d10Counts      = d10_cnt; 
      if(fCalculateCurrent){
	 // use your own calibration coefficients
	 ApplyCalibrationCoeff(pt);
      }else{
	 // use the values from the ROOT file
	 pt.u1Current    = u1_cur;
	 pt.unewCurrent  = unew_cur;
	 pt.d1Current    = d1_cur;
	 pt.d3Current    = d3_cur;
	 pt.d10Current   = d10_cur;
	 pt.dnewCurrent  = dnew_cur;
	 pt.unserCurrent = unser_cur;
      }
      if(arm.compare("Left")==0){
	 pt.time     = fLastTimeLeft + time; 
	 pt.time103kHz = fLastTimeLeft + time103kHz; 
	 pt.event    = fEvtCntrLeft;
	 if(runNumber==fLastRun){
	    pt.runEvent = fLastRunEvtCntLeft + i; 
	 }else{
	    pt.runEvent = i; 
	 }
	 fLeft.push_back(pt); 
	 fEvtCntrLeft++;
      }else if(arm.compare("sbs")==0){
	 pt.time       = fLastTimeSBS + time; 
	 pt.time103kHz = fLastTimeSBS + time103kHz; 
	 pt.event    = fEvtCntrSBS; 
	 if(runNumber==fLastRun){
	    pt.runEvent = fLastRunEvtCntSBS + i; 
	 }else{
	    pt.runEvent = i; 
	 }
	 fSBS.push_back(pt); 
	 fEvtCntrSBS++;
      }
   }

   // track number of events for this run
   if(arm.compare("sbs")==0)  fLastRunEvtCntSBS  = NN; 
   if(arm.compare("Left")==0) fLastRunEvtCntLeft = NN; 
 
   // track the last time registered 
   double lastTime=0; 
   if(arm.compare("Left")==0){
      if( IsBad(fLeft[NN-1].time) ){
         fLeft[NN-1].time = 0; 
	 std::cout << "[BCMManager::LoadDataFromTree]: WARNING! arm = " << arm << ", bad last time! Setting to zero." << std::endl;
      }
      fLastTimeLeft = fLeft[NN-1].time;
      lastTime = fLastTimeLeft; 
   }else if(arm.compare("SBS")==0){
      if( IsBad(fSBS[NN-1].time) ){
         fSBS[NN-1].time = 0; 
	 std::cout << "[BCMManager::LoadDataFromTree]: WARNING! arm = " << arm << ", bad last time! Setting to zero." << std::endl;
      } 
      fLastTimeSBS  = fSBS[NN-1].time;
      lastTime      = fLastTimeSBS; 
   }

   if(fIsDebug) std::cout << "The last time is: " << lastTime << std::endl;
 
   delete aTree; 
   delete ch;
   
   return 0; 
}
//______________________________________________________________________________
int BCMManager::LoadEPICSDataFromTree(const char *filePath,int runNumber){

   if(fIsDebug) std::cout << "[BCMManager::LoadEPICSDataFromTree]: Loading data from tree: E" << std::endl; 

   double epicsTime=0,IBC1H04CRCUR2=0,hac_bcm_average=0,halla_p=0;
   double IPM1H04A_XPOS=0,IPM1H04A_YPOS=0; 
   double IPM1H04E_XPOS=0,IPM1H04E_YPOS=0;
   double hac_bcm_dvm1_read=0,hac_bcm_dvm2_read=0; 
   double hac_bcm_dvm1_current=0,hac_bcm_dvm2_current=0; 

   TChain *ch = nullptr; 
   ch = new TChain("E"); 
   ch->Add(filePath);
   int NN = ch->GetEntries(); 

   if(ch==nullptr){
      return 1;
   }

   TTree *aTree = ch->GetTree();
   if(aTree==nullptr){
      return 1;
   }

   aTree->SetBranchAddress("halla_p"             ,&halla_p             ); 
   aTree->SetBranchAddress("IPM1H04A.XPOS"       ,&IPM1H04A_XPOS       );
   aTree->SetBranchAddress("IPM1H04A.YPOS"       ,&IPM1H04A_YPOS       );
   aTree->SetBranchAddress("IPM1H04E.XPOS"       ,&IPM1H04E_XPOS       );
   aTree->SetBranchAddress("IPM1H04E.YPOS"       ,&IPM1H04E_YPOS       );
   aTree->SetBranchAddress("hac_bcm_average"     ,&hac_bcm_average     );
   aTree->SetBranchAddress("hac_bcm_dvm1_read"   ,&hac_bcm_dvm1_read   );
   aTree->SetBranchAddress("hac_bcm_dvm2_read"   ,&hac_bcm_dvm2_read   );
   aTree->SetBranchAddress("hac_bcm_dvm1_current",&hac_bcm_dvm1_current);
   aTree->SetBranchAddress("hac_bcm_dvm2_current",&hac_bcm_dvm1_current);
   aTree->SetBranchAddress("IBC1H04CRCUR2"       ,&IBC1H04CRCUR2       );
   aTree->SetBranchAddress("timestamp"           ,&epicsTime           );

   epicsData_t pt; 
   for(int i=0;i<NN;i++){
      aTree->GetEntry(i); 
      pt.event           = fEvtCntrEPICS; 
      pt.runNumber       = runNumber;
      pt.time            = epicsTime;
      if(runNumber==fLastRun){
	 pt.runEvent = fLastRunEvtCntEPICS + i; 
      }else{
	 pt.runEvent = i; 
      }
      pt.halla_p              = halla_p;
      pt.hac_bcm_average      = hac_bcm_average; 
      pt.hac_bcm_dvm1_read    = hac_bcm_dvm1_read; 
      pt.hac_bcm_dvm2_read    = hac_bcm_dvm2_read; 
      pt.hac_bcm_dvm1_current = hac_bcm_dvm1_current; 
      pt.hac_bcm_dvm2_current = hac_bcm_dvm2_current; 
      pt.IBC1H04CRCUR2        = IBC1H04CRCUR2; 
      pt.IPM1H04A_XPOS        = IPM1H04A_XPOS;   
      pt.IPM1H04A_YPOS        = IPM1H04A_YPOS;   
      pt.IPM1H04E_XPOS        = IPM1H04E_XPOS;   
      pt.IPM1H04E_YPOS        = IPM1H04E_YPOS; 
      fEPICS.push_back(pt);  
      fEvtCntrEPICS++; 
   } 

   // track number of events for this run
   fLastRunEvtCntEPICS = NN; 

   delete aTree;
   delete ch; 

   return 0;
}
//______________________________________________________________________________
int BCMManager::LoadCalibrationCoefficients(const char *dirName){
   fCalculateCurrent = true;
   // load calibration coefficients from a specific directory/production set 
   const int N = 7; 
   std::string bcm[N] = {"unser","u1","unew","d1","d3","d10","dnew"};

   char inpath[200]; 

   for(int i=0;i<N;i++){
      sprintf(inpath,"%s/calib-coeff_%s.csv",dirName,bcm[i].c_str()); 
      LoadCalibrationCoefficients(bcm[i].c_str(),inpath); 
   } 
   return 0;
}
//______________________________________________________________________________
int BCMManager::LoadCalibrationCoefficients(const char *type,const char *filePath){
   std::string Type = type; 

   CSVManager *csv = new CSVManager();
   csv->ReadFile(filePath,true); 

   std::vector<std::string> dev;
   std::vector<int> runMin,runMax;  
   std::vector<double> ped,pedErr,offset,offsetErr,gain,gainErr;

   csv->GetColumn_byName_str("dev",dev);
   csv->GetColumn_byName<int>("runMin",runMin);  
   csv->GetColumn_byName<int>("runMax",runMax);  
   csv->GetColumn_byName<double>("pedestal"   ,ped); 
   csv->GetColumn_byName<double>("pedestalErr",pedErr); 
   csv->GetColumn_byName<double>("offset"     ,offset); 
   csv->GetColumn_byName<double>("offsetErr"  ,offsetErr); 
   csv->GetColumn_byName<double>("gain"       ,gain); 
   csv->GetColumn_byName<double>("gainErr"    ,gainErr); 

   calibCoeff_t cc; 
   const int ND = dev.size();
   for(int i=0;i<ND;i++){
      cc.dev         = dev[i];
      cc.runMin      = runMin[i]; 
      cc.runMax      = runMax[i]; 
      cc.pedestal    = ped[i];  
      cc.pedestalErr = pedErr[i];  
      cc.offset      = offset[i];  
      cc.offsetErr   = offsetErr[i];  
      cc.slope       = gain[i];  
      cc.slopeErr    = gainErr[i];  
      if( Type.compare("unser")==0) fccUnser.push_back(cc); 
      if( Type.compare("u1")==0   ) fccU1.push_back(cc); 
      if( Type.compare("unew")==0 ) fccUnew.push_back(cc); 
      if( Type.compare("d1")==0   ) fccD1.push_back(cc); 
      if( Type.compare("d3")==0   ) fccD3.push_back(cc); 
      if( Type.compare("d10")==0  ) fccD10.push_back(cc); 
      if( Type.compare("dnew")==0 ) fccDnew.push_back(cc); 
   }
   delete csv;
   return 0; 
}
//______________________________________________________________________________
int BCMManager::ApplyCalibrationCoeff(scalerData_t &data){
   // apply calibration coefficients to all data
   std::vector<std::string> dev;
   dev.push_back("unser"); 
   dev.push_back("u1"); 
   dev.push_back("unew"); 
   dev.push_back("d1"); 
   dev.push_back("d3"); 
   dev.push_back("d10"); 
   dev.push_back("dnew");

   double PED=0,OFFSET=0,GAIN=0;
   int rc=0;  

   calibCoeff_t ccData;
 
   const int ND = dev.size();
   for(int i=0;i<ND;i++){
      rc = GetCalibrationCoeff(data.runNumber,dev[i],ccData);
      PED = ccData.pedestal; OFFSET = ccData.offset; GAIN = ccData.slope;
      if(rc==0){  
         // if rc = 0, we found the calibration coefficients; apply them
	 if(dev[i].compare("unser")==0){
	    data.unserCurrent = (data.unserRate - PED - OFFSET)/GAIN;
            data.unserRate_ps = data.unserRate - PED; 
	 }else if(dev[i].compare("u1")==0){
	    data.u1Current = (data.u1Rate - PED - OFFSET)/GAIN;
            data.u1Rate_ps = data.u1Rate - PED; 
	 }else if(dev[i].compare("unew")==0){
	    data.unewCurrent = (data.unewRate - PED - OFFSET)/GAIN;
            data.unewRate_ps = data.unewRate - PED; 
	 }else if(dev[i].compare("d1")==0){
	    data.d1Current = (data.d1Rate - PED - OFFSET)/GAIN;
            data.d1Rate_ps = data.d1Rate - PED; 
	 }else if(dev[i].compare("d3")==0){
	    data.d3Current = (data.d3Rate - PED - OFFSET)/GAIN;
            data.d3Rate_ps = data.d3Rate - PED; 
	 }else if(dev[i].compare("d10")==0){
	    data.d10Current = (data.d10Rate - PED - OFFSET)/GAIN;
            data.d10Rate_ps = data.d10Rate - PED; 
	 }else if(dev[i].compare("dnew")==0){
	    data.dnewCurrent = (data.dnewRate - PED - OFFSET)/GAIN;
            data.dnewRate_ps = data.dnewRate - PED; 
	 }
      }
   } 
   return 0;
}
//______________________________________________________________________________
int BCMManager::GetCalibrationCoeff(int run,std::string dev,calibCoeff_t &data){
   // get calibration coefficients by device, and based on run number
   int NGRP=0;
   std::vector<calibCoeff_t> cc; 
   if(dev.compare("unser")==0){
      NGRP = fccUnser.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccUnser[i]); 
   }else if(dev.compare("u1")==0){
      NGRP = fccU1.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccU1[i]); 
   }else if(dev.compare("unew")==0){
      NGRP = fccUnew.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccUnew[i]); 
   }else if(dev.compare("d1")==0){
      NGRP = fccD1.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccD1[i]); 
   }else if(dev.compare("d3")==0){
      NGRP = fccD3.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccD3[i]); 
   }else if(dev.compare("d10")==0){
      NGRP = fccD10.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccD10[i]); 
   }else if(dev.compare("dnew")==0){
      NGRP = fccDnew.size();
      for(int i=0;i<NGRP;i++) cc.push_back(fccDnew[i]); 
   }
   // find coeffs by run group 
   for(int i=0;i<NGRP;i++){
      if((run>=cc[i].runMin)&&(run<=cc[i].runMax)){
	 // matched run! this is what we need 
	 data = cc[i];  
      }
   }
   // now check the gain to avoid div by zero error 
   if(data.slope==0){
      if(fVerbosity>5) std::cout << "[BCMManager::GetCalibrationCoeff]: WARNING! Invalid gain = 0! Defaulting to gain = 1" << std::endl;
      data.slope = 1; 
   }
   return 0;
}
//______________________________________________________________________________
bool BCMManager::IsBad(double v){
   return std::isinf(v) || std::isnan(v); 
}
//______________________________________________________________________________
TH1F * BCMManager::GetTH1F(const char *arm,const char *var_name,int NBin,double min,double max){
   std::vector<double> x;
   GetVector(arm,var_name,x); 

   TString name = Form("%s.%s",arm,var_name);
   TH1F *h = new TH1F(name,name,NBin,min,max); 
   const int N = x.size();
   for(int i=0;i<N;i++) h->Fill(x[i]); 

   return h; 
}
//______________________________________________________________________________
TH2F * BCMManager::GetTH2F(const char *arm,const char *xAxis,const char *yAxis,
                           int NXBin,double xMin,double xMax,int NYBin,double yMin,double yMax){

   std::vector<double> x,y;
   GetVector(arm,xAxis,x); 
   GetVector(arm,yAxis,y); 

   TString name = Form("%s.%s:%s.%s",arm,yAxis,arm,xAxis);
   TH2F *h = new TH2F(name,name,NXBin,xMin,xMax,NYBin,yMin,yMax); 
   const int N = x.size();
   for(int i=0;i<N;i++) h->Fill(x[i],y[i]); 

   return h;
}
//______________________________________________________________________________
TGraph * BCMManager::GetTGraph(const char *arm,const char *xAxis,const char *yAxis){
   // get a TGraph object 
   std::string ARM = arm; 
   TString xAxisName = Form("%s",xAxis);
   TString yAxisName = Form("%s",yAxis); 

   if(fIsDebug) std::cout << xAxisName << " " << yAxisName << std::endl;

   std::vector<double> x,y; 
   GetVector(arm,xAxisName.Data(),x); 
   GetVector(arm,yAxisName.Data(),y);

   const int N = x.size(); 
   TGraph *g = new TGraph(N,&x[0],&y[0]); 
   return g; 
}
//______________________________________________________________________________
int BCMManager::GetVector_scaler(const char *arm,int run,std::vector<scalerData_t> &data){
   // return a vector of the scalerData for a given run for a given run for a given run for a given run 
   std::string armName = arm; 
 
   int NN=0;
   if(armName.compare("Left")==0) NN = fLeft.size(); 
   if(armName.compare("sbs")==0)  NN = fSBS.size(); 

   for(int i=0;i<NN;i++){
      if(armName.compare("Left")==0){
	 if(run==fLeft[i].runNumber){
	    data.push_back(fLeft[i]);
         } 
      }else if(armName.compare("sbs")==0){
	 if(run==fSBS[i].runNumber){
	    data.push_back(fSBS[i]);
         } 
      }
   } 
   return 0;
}
//______________________________________________________________________________
int BCMManager::GetVector_scaler(const char *arm,std::vector<scalerData_t> &data){
   // return a vector of the scalerData 
   std::string armName = arm; 
 
   int NN=0;
   if(armName.compare("Left")==0) NN = fLeft.size(); 
   if(armName.compare("sbs")==0)  NN = fSBS.size(); 

   for(int i=0;i<NN;i++){
      if(armName.compare("Left")==0){
	 data.push_back(fLeft[i]); 
      }else if(armName.compare("sbs")==0){
	 data.push_back(fSBS[i]); 
      }
   } 
   return 0;
}
//______________________________________________________________________________
int BCMManager::GetVector_epics(std::vector<epicsData_t> &data){
   // return a vector of the epicsData 
   int NN = fEPICS.size();
   for(int i=0;i<NN;i++){
      data.push_back(fEPICS[i]); 
   } 
   return 0;
}
//______________________________________________________________________________
int BCMManager::GetVector_epics(int run,std::vector<epicsData_t> &data){
   // return a vector of the epicsData 
   int NN = fEPICS.size();
   for(int i=0;i<NN;i++){
      if(run==fEPICS[i].runNumber) data.push_back(fEPICS[i]); 
   } 
   return 0;
}
//______________________________________________________________________________
int BCMManager::GetVector(const char *arm,const char *var,std::vector<double> &v){
   // fill a vector with the variable 

   std::string armName = arm;
   std::string varName = var;

   // if(fIsDebug) std::cout << "arm = " << arm << ", var = " << var << std::endl;

   int NN=0;
   if(armName.compare("Left")==0) NN = fLeft.size(); 
   if(armName.compare("sbs")==0)  NN = fSBS.size(); 
   if(armName.compare("E")==0)    NN = fEPICS.size(); 

   double val=0;
   for(int i=0;i<NN;i++){
      val = 0;
      if(armName.compare("Left")==0) val = fLeft[i].getValue(varName); 
      if(armName.compare("sbs")==0)  val = fSBS[i].getValue(varName); 
      if(armName.compare("E")==0)    val = fEPICS[i].getValue(varName); 
      v.push_back(val);
   }
  
   int NV = v.size();
   if(NV==0) std::cout << "[BCMManager::GetVector]: No events for arm = " << armName << ", variable = " << varName << "!" << std::endl;

   return 0;
}
//______________________________________________________________________________
void BCMManager::Print(const char *arm){
   // print data to screen
   std::string ARM = arm; 
   const int NL = fLeft.size();
   const int NS = fSBS.size();
   const int NE = fEPICS.size();
   int N=0;
   if(ARM.compare("Left")==0) N = NL; 
   if(ARM.compare("sbs")==0)  N = NS; 
   if(ARM.compare("E")==0)    N = NE; 

   for(int i=0;i<N;i++){
      if(ARM.compare("Left")==0){
	 fLeft[i].Print("rate"); 
      }else if(ARM.compare("sbs")==0){
	 fSBS[i].Print("rate"); 
      }else if(ARM.compare("E")==0){
	 fEPICS[i].Print();
      } 
   }
}

#include "../include/BeamManager.h"
//______________________________________________________________________________
BeamManager::BeamManager(const char *filePath,const char *ccDirPath,bool isDebug){
   fIsDebug          = isDebug; 
   fCalculateCurrent = false; 
   fVerbosity        = 0;  
   fEvtCntrLeft      = 0;  
   fEvtCntrSBS       = 0;  
   fEvtCntrEPICS     = 0; 
   fLastTimeLeft     = 0; 
   fLastTimeSBS      = 0;

   fDeltaS_bpm       = 4.083; // distance between BPMA and BPMB in meters 
   fS_bpmA           = 5.969; // distance from BPMA to target pivot         

   std::string path = filePath; 
   if(path.compare("NONE")!=0){
      LoadFile(filePath);
   }

   // std::string cc_dirPath = ccDirPath; 
   // if(cc_dirPath.compare("NONE")!=0){
   //    LoadCalibrationCoefficients(ccDirPath);
   // } 
 
}
//______________________________________________________________________________
BeamManager::~BeamManager(){
   Clear();
}
//______________________________________________________________________________
void BeamManager::Clear(){
   fLeft.clear();
   fSBS.clear();
   fEPICS.clear();
   fRunList.clear(); 
   fEvtCntrLeft  = 0;
   fEvtCntrSBS   = 0;
   fEvtCntrEPICS = 0;
   fLastTimeLeft = 0; 
   fLastTimeSBS  = 0; 
}
//______________________________________________________________________________
int BeamManager::LoadFile(const char *filePath,int runNumber){
   // Attach trees for LHRS and SBS
   int rc = CheckFile(filePath); 
   if(rc!=0){
      std::cout << "[BeamManager::LoadFile]: Cannot open the file for run " << runNumber << std::endl;
      return rc;
   }
   // file is valid, read in the data from the trees 
   int rc_lhrs  = LoadDataFromTree(filePath,"Lrb"  ,runNumber); 
   int rc_sbs   = LoadDataFromTree(filePath,"SBSrb",runNumber); 
   int rc_epics = LoadEPICSDataFromTree(filePath,runNumber); 

   // add the run to the run list
   // must at least have the SBS branch to be considered good 
   if(rc_sbs==0){
      std::cout << "[BeamManager::LoadFile]: Loaded run " << runNumber << std::endl;
      fRunList.push_back(runNumber);
   }

   return rc; 
}
//______________________________________________________________________________
int BeamManager::CheckFile(const char *filePath){
   TFile *myFile = new TFile(filePath);

   Long64_t bytesRead = myFile->GetBytesRead();

   int rc=1; 
   if(bytesRead!=0){
      // file has non-zero size
      rc = 0;
   } 
   return rc;
}
//______________________________________________________________________________
int BeamManager::LoadDataFromTree(const char *filePath,const char *branchName,int runNumber){

   if(fIsDebug) std::cout << "[BeamManager::LoadDataFromTree]: Loading data from tree: " << branchName << std::endl;

   double raster1_rawcur_x=0,raster1_rawcur_y=0; 
   double raster2_rawcur_x=0,raster2_rawcur_y=0;
   double bpmA_rawcur_1=0,bpmA_rawcur_2=0,bpmA_rawcur_3=0,bpmA_rawcur_4=0; 
   double bpmB_rawcur_1=0,bpmB_rawcur_2=0,bpmB_rawcur_3=0,bpmB_rawcur_4=0; 
   double bpmA_x=0,bpmA_y=0,bpmB_x=0,bpmB_y=0;  

   TChain *ch = nullptr; 
   ch = new TChain("T"); 
   ch->Add(filePath);
   int NN = ch->GetEntries();

   if(ch==nullptr){
      return 1;
   } 
 
   TTree *aTree = ch->GetTree();
   if(aTree==nullptr) return 1; 

   aTree->SetBranchAddress(Form("%s.Raster.rawcur.x" ,branchName),&raster1_rawcur_x );
   aTree->SetBranchAddress(Form("%s.Raster.rawcur.y" ,branchName),&raster1_rawcur_y );
   aTree->SetBranchAddress(Form("%s.Raster2.rawcur.x",branchName),&raster2_rawcur_x );
   aTree->SetBranchAddress(Form("%s.Raster2.rawcur.y",branchName),&raster2_rawcur_y );
   aTree->SetBranchAddress(Form("%s.BPMA.rawcur.1"   ,branchName),&bpmA_rawcur_1 );
   aTree->SetBranchAddress(Form("%s.BPMA.rawcur.2"   ,branchName),&bpmA_rawcur_2 );
   aTree->SetBranchAddress(Form("%s.BPMA.rawcur.3"   ,branchName),&bpmA_rawcur_3 );
   aTree->SetBranchAddress(Form("%s.BPMA.rawcur.4"   ,branchName),&bpmA_rawcur_4 );
   aTree->SetBranchAddress(Form("%s.BPMB.rawcur.1"   ,branchName),&bpmB_rawcur_1 );
   aTree->SetBranchAddress(Form("%s.BPMB.rawcur.2"   ,branchName),&bpmB_rawcur_2 );
   aTree->SetBranchAddress(Form("%s.BPMB.rawcur.3"   ,branchName),&bpmB_rawcur_3 );
   aTree->SetBranchAddress(Form("%s.BPMB.rawcur.4"   ,branchName),&bpmB_rawcur_4 );
   aTree->SetBranchAddress(Form("%s.BPMA.x"          ,branchName),&bpmA_x );
   aTree->SetBranchAddress(Form("%s.BPMA.y"          ,branchName),&bpmA_y );
   aTree->SetBranchAddress(Form("%s.BPMB.x"          ,branchName),&bpmB_x );
   aTree->SetBranchAddress(Form("%s.BPMB.y"          ,branchName),&bpmB_y );

   beamData_t pt; 

   std::string Branch = branchName; 

   for(int i=0;i<NN;i++){
      aTree->GetEntry(i);
      // if(time_den!=0){
      //    time = time_num/time_den; 
      // }else{
      //    time = 0;  
      // } 
      // pt.time_num    = time_num;  
      // pt.time_den    = time_den; 
      pt.arm           = branchName; 
      pt.runNumber     = runNumber;
      pt.bpmA_x        = bpmA_x; 
      pt.bpmA_y        = bpmA_y; 
      pt.bpmB_x        = bpmB_x; 
      pt.bpmB_y        = bpmB_y; 
      pt.bpmA_rawcur_1 = bpmA_rawcur_1; 
      pt.bpmA_rawcur_2 = bpmA_rawcur_2; 
      pt.bpmA_rawcur_3 = bpmA_rawcur_3; 
      pt.bpmA_rawcur_4 = bpmA_rawcur_4; 
      pt.bpmB_rawcur_1 = bpmB_rawcur_1; 
      pt.bpmB_rawcur_2 = bpmB_rawcur_2; 
      pt.bpmB_rawcur_3 = bpmB_rawcur_3; 
      pt.bpmB_rawcur_4 = bpmB_rawcur_4; 
      pt.raster1_rawcur_x = raster1_rawcur_x; 
      pt.raster1_rawcur_y = raster1_rawcur_y; 
      pt.raster2_rawcur_x = raster2_rawcur_x; 
      pt.raster2_rawcur_y = raster2_rawcur_y;
      // compute beam position at the target [in mm]  
      pt.target_x = 1000.*( bpmA_x + (bpmB_x-bpmA_x)*(fS_bpmA/fDeltaS_bpm) );
      pt.target_y = 1000.*( bpmA_y + (bpmB_y-bpmA_y)*(fS_bpmA/fDeltaS_bpm) );
      // if(fCalculateCurrent){
      //    ApplyCalibrationCoeff(pt);
      // }
      if(Branch.compare("Lrb")==0){
	 // pt.time  = fLastTimeLeft + time; 
	 pt.event    = fEvtCntrLeft; 
	 pt.runEvent = i;  
	 fLeft.push_back(pt); 
	 fEvtCntrLeft++;
      }else if(Branch.compare("SBSrb")==0){
	 // pt.time  = fLastTimeSBS + time; 
	 pt.event    = fEvtCntrSBS; 
	 pt.runEvent = i;  
	 fSBS.push_back(pt); 
	 fEvtCntrSBS++;
      }
   }
  
   // track the last time registered 
   // double lastTime=0; 
   // if(Branch.compare("Lrb")==0){
   //    if( IsBad(fLeft[NN-1].time) ){
   //       fLeft[NN-1].time = 0; 
   //       std::cout << "[BeamManager::LoadDataFromTree]: WARNING! arm = " << branchName << ", bad last time! Setting to zero." << std::endl;
   //    }
   //    fLastTimeLeft = fLeft[NN-1].time;
   //    lastTime = fLastTimeLeft; 
   // }else if(arm.compare("SBSrb")==0){
   //    if( IsBad(fSBS[NN-1].time) ){
   //       std::cout << "[BeamManager::LoadDataFromTree]: WARNING! arm = " << branchName << ", bad last time! Setting to zero." << std::endl;
   //       fSBS[NN-1].time = 0; 
   //    } 
   //    fLastTimeSBS  = fSBS[NN-1].time;
   //    lastTime      = fLastTimeSBS; 
   // }
   // if(fIsDebug) std::cout << "The last time is: " << lastTime << std::endl;
 
   delete aTree; 
   delete ch;
   
   return 0; 
}
//______________________________________________________________________________
int BeamManager::LoadEPICSDataFromTree(const char *filePath,int runNumber){

   if(fIsDebug) std::cout << "[BeamManager::LoadEPICSDataFromTree]: Loading data from tree: E" << std::endl; 

   double epicsTime=0,IBC1H04CRCUR2=0,hac_bcm_average=0;
   double IPM1H04A_XPOS=0,IPM1H04A_YPOS=0; 
   double IPM1H04E_XPOS=0,IPM1H04E_YPOS=0; 

   TChain *ch = nullptr; 
   ch = new TChain("E"); 
   ch->Add(filePath);

   if(ch==nullptr){
      return 1;
   }

   int NL = ch->GetEntries(); 
   TTree *aTree = ch->GetTree();
   if(aTree==nullptr){
      return 1;
   }

   aTree->SetBranchAddress("hac_bcm_average",&hac_bcm_average);
   aTree->SetBranchAddress("IPM1H04A.XPOS"  ,&IPM1H04A_XPOS  );
   aTree->SetBranchAddress("IPM1H04A.YPOS"  ,&IPM1H04A_YPOS  );
   aTree->SetBranchAddress("IPM1H04E.XPOS"  ,&IPM1H04E_XPOS  );
   aTree->SetBranchAddress("IPM1H04E.YPOS"  ,&IPM1H04E_YPOS  );
   aTree->SetBranchAddress("IBC1H04CRCUR2"  ,&IBC1H04CRCUR2  );
   aTree->SetBranchAddress("timestamp"      ,&epicsTime      );

   epicsData_t pt; 
   for(int i=0;i<NL;i++){
      aTree->GetEntry(i); 
      pt.event           = fEvtCntrEPICS; 
      pt.runNumber       = runNumber;
      pt.time            = epicsTime;
      pt.hac_bcm_average = hac_bcm_average; 
      pt.IBC1H04CRCUR2   = IBC1H04CRCUR2; 
      pt.IPM1H04A_XPOS   = IPM1H04A_XPOS;   
      pt.IPM1H04A_YPOS   = IPM1H04A_YPOS;   
      pt.IPM1H04E_XPOS   = IPM1H04E_XPOS;   
      pt.IPM1H04E_YPOS   = IPM1H04E_YPOS; 
      fEPICS.push_back(pt);  
      fEvtCntrEPICS++; 
   }

   delete aTree;
   delete ch; 

   return 0;
}
// //______________________________________________________________________________
// int BeamManager::LoadCalibrationCoefficients(const char *dirName){
//    fCalculateCurrent = true;
//    // load calibration coefficients from a specific directory/production set 
//    const int N = 7; 
//    std::string bcm[N] = {"unser","u1","unew","d1","d3","d10","dnew"};
// 
//    char inpath[200]; 
// 
//    for(int i=0;i<N;i++){
//       sprintf(inpath,"%s/calib-coeff_%s.csv",dirName,bcm[i].c_str()); 
//       LoadCalibrationCoefficients(bcm[i].c_str(),inpath); 
//    } 
//    return 0;
// }
// //______________________________________________________________________________
// int BeamManager::LoadCalibrationCoefficients(const char *type,const char *filePath){
//    std::string Type = type; 
// 
//    CSVManager *csv = new CSVManager();
//    csv->ReadFile(filePath,true); 
// 
//    std::vector<std::string> dev;
//    std::vector<int> runMin,runMax;  
//    std::vector<double> ped,pedErr,offset,offsetErr,gain,gainErr;
// 
//    csv->GetColumn_byName_str("dev",dev);
//    csv->GetColumn_byName<int>("runMin",runMin);  
//    csv->GetColumn_byName<int>("runMax",runMax);  
//    csv->GetColumn_byName<double>("pedestal"   ,ped); 
//    csv->GetColumn_byName<double>("pedestalErr",pedErr); 
//    csv->GetColumn_byName<double>("offset"     ,offset); 
//    csv->GetColumn_byName<double>("offsetErr"  ,offsetErr); 
//    csv->GetColumn_byName<double>("gain"       ,gain); 
//    csv->GetColumn_byName<double>("gainErr"    ,gainErr); 
// 
//    calibCoeff_t cc; 
//    const int ND = dev.size();
//    for(int i=0;i<ND;i++){
//       cc.dev         = dev[i];
//       cc.runMin      = runMin[i]; 
//       cc.runMax      = runMax[i]; 
//       cc.pedestal    = ped[i];  
//       cc.pedestalErr = pedErr[i];  
//       cc.offset      = offset[i];  
//       cc.offsetErr   = offsetErr[i];  
//       cc.slope       = gain[i];  
//       cc.slopeErr    = gainErr[i];  
//       if( Type.compare("unser")==0) fccUnser.push_back(cc); 
//       if( Type.compare("u1")==0   ) fccU1.push_back(cc); 
//       if( Type.compare("unew")==0 ) fccUnew.push_back(cc); 
//       if( Type.compare("d1")==0   ) fccD1.push_back(cc); 
//       if( Type.compare("d3")==0   ) fccD3.push_back(cc); 
//       if( Type.compare("d10")==0  ) fccD10.push_back(cc); 
//       if( Type.compare("dnew")==0 ) fccDnew.push_back(cc); 
//    }
//    delete csv;
//    return 0; 
// }
// //______________________________________________________________________________
// int BeamManager::ApplyCalibrationCoeff(beamData_t &data){
//    // apply calibration coefficients to all data
//    std::vector<std::string> dev;
//    dev.push_back("unser"); 
//    dev.push_back("u1"); 
//    dev.push_back("unew"); 
//    dev.push_back("d1"); 
//    dev.push_back("d3"); 
//    dev.push_back("d10"); 
//    dev.push_back("dnew");
// 
//    double PED=0,OFFSET=0,GAIN=0;
//    int rc=0;  
// 
//    calibCoeff_t ccData;
//  
//    const int ND = dev.size();
//    for(int i=0;i<ND;i++){
//       rc = GetCalibrationCoeff(data.runNumber,dev[i],ccData);
//       PED = ccData.pedestal; OFFSET = ccData.offset; GAIN = ccData.slope;
//       if(rc==0){  
//          // if rc = 0, we found the calibration coefficients; apply them
// 	 if(dev[i].compare("unser")==0){
// 	    data.unserCurrent = (data.unserRate - PED - OFFSET)/GAIN;
// 	 }else if(dev[i].compare("u1")==0){
// 	    data.u1Current = (data.u1Rate - PED - OFFSET)/GAIN;
// 	 }else if(dev[i].compare("unew")==0){
// 	    data.unewCurrent = (data.unewRate - PED - OFFSET)/GAIN;
// 	 }else if(dev[i].compare("d1")==0){
// 	    data.d1Current = (data.d1Rate - PED - OFFSET)/GAIN;
// 	 }else if(dev[i].compare("d3")==0){
// 	    data.d3Current = (data.d3Rate - PED - OFFSET)/GAIN;
// 	 }else if(dev[i].compare("d10")==0){
// 	    data.d10Current = (data.d10Rate - PED - OFFSET)/GAIN;
// 	 }else if(dev[i].compare("dnew")==0){
// 	    data.dnewCurrent = (data.dnewRate - PED - OFFSET)/GAIN;
// 	 }
//       }
//    } 
//    return 0;
// }
// //______________________________________________________________________________
// int BeamManager::GetCalibrationCoeff(int run,std::string dev,calibCoeff_t &data){
//    // get calibration coefficients by device, and based on run number
//    int NGRP=0;
//    std::vector<calibCoeff_t> cc; 
//    if(dev.compare("unser")==0){
//       NGRP = fccUnser.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccUnser[i]); 
//    }else if(dev.compare("u1")==0){
//       NGRP = fccU1.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccU1[i]); 
//    }else if(dev.compare("unew")==0){
//       NGRP = fccUnew.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccUnew[i]); 
//    }else if(dev.compare("d1")==0){
//       NGRP = fccD1.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccD1[i]); 
//    }else if(dev.compare("d3")==0){
//       NGRP = fccD3.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccD3[i]); 
//    }else if(dev.compare("d10")==0){
//       NGRP = fccD10.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccD10[i]); 
//    }else if(dev.compare("dnew")==0){
//       NGRP = fccDnew.size();
//       for(int i=0;i<NGRP;i++) cc.push_back(fccDnew[i]); 
//    }
//    // find coeffs by run group 
//    for(int i=0;i<NGRP;i++){
//       if((run>=cc[i].runMin)&&(run<=cc[i].runMax)){
// 	 // matched run! this is what we need 
// 	 data = cc[i];  
//       }
//    }
//    // now check the gain to avoid div by zero error 
//    if(data.slope==0){
//       if(fVerbosity>5) std::cout << "[BeamManager::GetCalibrationCoeff]: WARNING! Invalid gain = 0! Defaulting to gain = 1" << std::endl;
//       data.slope = 1; 
//    }
//    return 0;
// }
//______________________________________________________________________________
bool BeamManager::IsBad(double v){
   return std::isinf(v) || std::isnan(v); 
}
//______________________________________________________________________________
TH1F * BeamManager::GetTH1F(const char *arm,const char *var_name,int NBin,double min,double max){
   std::vector<double> x;
   GetVector(arm,var_name,x); 

   TString name = Form("%s.%s",arm,var_name);
   TH1F *h = new TH1F(name,name,NBin,min,max); 
   const int N = x.size();
   for(int i=0;i<N;i++) h->Fill(x[i]); 

   return h; 
}
//______________________________________________________________________________
TH2F * BeamManager::GetTH2F(const char *arm,const char *xAxis,const char *yAxis,
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
TGraph * BeamManager::GetTGraph(const char *arm,const char *xAxis,const char *yAxis){
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
int BeamManager::GetVector_beam(const char *arm,std::vector<beamData_t> &data){
   // return a vector of the beamData 
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
int BeamManager::GetVector_epics(std::vector<epicsData_t> &data){
   // return a vector of the epicsData 
   int NN = fEPICS.size();
   for(int i=0;i<NN;i++){
      data.push_back(fEPICS[i]); 
   } 
   return 0;
}
//______________________________________________________________________________
int BeamManager::GetVector(const char *arm,const char *var,std::vector<double> &v){
   // fill a vector with the variable 

   std::string armName = arm;
   std::string varName = var;

   // if(fIsDebug) std::cout << "arm = " << arm << ", var = " << var << std::endl;

   int NN=0;
   if(armName.compare("Lrb")==0)    NN = fLeft.size(); 
   if(armName.compare("SBSrb")==0)  NN = fSBS.size(); 
   if(armName.compare("E")==0)      NN = fEPICS.size(); 

   double val=0;
   for(int i=0;i<NN;i++){
      val = 0;
      if(armName.compare("Lrb")==0)   val = fLeft[i].getValue(varName); 
      if(armName.compare("SBSrb")==0) val = fSBS[i].getValue(varName); 
      if(armName.compare("E")==0)     val = fEPICS[i].getValue(varName); 
      v.push_back(val);
   }
  
   int NV = v.size();
   if(NV==0) std::cout << "[BeamManager::GetVector]: No events for arm = " << armName << ", variable = " << varName << "!" << std::endl;

   return 0;
}
//______________________________________________________________________________
void BeamManager::Print(const char *arm){
   // print data to screen
   std::string ARM = arm; 
   const int NL = fLeft.size();
   const int NS = fSBS.size();
   const int NE = fEPICS.size();
   int N=0;
   if(ARM.compare("Lrb")==0)   N = NL; 
   if(ARM.compare("SBSrb")==0) N = NS; 
   if(ARM.compare("E")==0)     N = NE; 

   for(int i=0;i<N;i++){
      if(ARM.compare("Lrb")==0){
	 fLeft[i].Print(); 
      }else if(ARM.compare("SBSrb")==0){
	 fSBS[i].Print(); 
      }else if(ARM.compare("E")==0){
	 fEPICS[i].Print();
      } 
   }
}

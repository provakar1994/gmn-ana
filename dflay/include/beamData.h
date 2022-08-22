#ifndef BEAM_DATA_H
#define BEAM_DATA_H

// a struct to organize beam data (BPMs, raster, helicity)

typedef struct beamData {
   std::string arm;               // Lrb or SBSrb 
   std::string info;              // user notes
   unsigned int timestamp;        // UTC timestamp 
   double time;                   // time 
   double raster1_rawcur_x;       // raster raw current in x (coil for By) (upstream)   
   double raster1_rawcur_y;       // raster raw current in y (coil for Bx) (upstream)
   double raster2_rawcur_x;       // raster raw current in x (coil for By) (downstream) 
   double raster2_rawcur_y;       // raster raw current in x (coil for Bx) (downstream)
   double bpmA_rawcur_1;          // beam position current 1, x+ (upstream)       
   double bpmA_rawcur_2;          // beam position current 2, x- (upstream)         
   double bpmA_rawcur_3;          // beam position current 3, y+ (upstream)         
   double bpmA_rawcur_4;          // beam position current 4, y- (upstream)         
   double bpmB_rawcur_1;          // beam position current 1, x+ (downstream)        
   double bpmB_rawcur_2;          // beam position current 2, x- (downstream)   
   double bpmB_rawcur_3;          // beam position current 3, y+ (downstream)   
   double bpmB_rawcur_4;          // beam position current 4, y- (downstream)   
   double bpmA_x;                 // beam position in x (upstream)       
   double bpmA_y;                 // beam position in y (upstream)         
   double bpmB_x;                 // beam position in x (downstream)        
   double bpmB_y;                 // beam position in y (downstream) 
   double target_x;               // beam position at target in x (using BPMs)  
   double target_y;               // beam position at target in y (using BPMs)  
   int runNumber;
   int event; 
   int runEvent;
   // constructor 
   beamData(): 
      arm("NONE"),info("NONE"),
      timestamp(0),time(0),
      raster1_rawcur_x(0),raster1_rawcur_y(0),raster2_rawcur_x(0),raster2_rawcur_y(0),
      bpmA_rawcur_1(0),bpmA_rawcur_2(0),bpmA_rawcur_3(0),bpmA_rawcur_4(0), 
      bpmB_rawcur_1(0),bpmB_rawcur_2(0),bpmB_rawcur_3(0),bpmB_rawcur_4(0), 
      bpmA_x(0),bpmA_y(0),bpmB_x(0),bpmB_y(0),
      target_x(0),target_y(0),
      runNumber(0),event(0),runEvent(0) 
   {}  
   double getValue(std::string name){
      double val=0;
      if(name.compare("event")==0)            val = (double)event; 
      if(name.compare("runEvent")==0)         val = (double)runEvent; 
      if(name.compare("run")==0)              val = (double)runNumber; 
      if(name.compare("time")==0)             val = (double)timestamp; 
      if(name.compare("timestamp")==0)        val = time; 
      if(name.compare("Raster1.rawcur.x")==0) val = raster1_rawcur_x; 
      if(name.compare("Raster1.rawcur.y")==0) val = raster1_rawcur_y; 
      if(name.compare("Raster2.rawcur.x")==0) val = raster2_rawcur_x; 
      if(name.compare("Raster2.rawcur.y")==0) val = raster2_rawcur_y; 
      if(name.compare("BPMA.rawcur.1")==0)    val = bpmA_rawcur_1; 
      if(name.compare("BPMA.rawcur.2")==0)    val = bpmA_rawcur_2; 
      if(name.compare("BPMA.rawcur.3")==0)    val = bpmA_rawcur_3; 
      if(name.compare("BPMA.rawcur.4")==0)    val = bpmA_rawcur_4; 
      if(name.compare("BPMB.rawcur.1")==0)    val = bpmB_rawcur_1; 
      if(name.compare("BPMB.rawcur.2")==0)    val = bpmB_rawcur_2; 
      if(name.compare("BPMB.rawcur.3")==0)    val = bpmB_rawcur_3; 
      if(name.compare("BPMB.rawcur.4")==0)    val = bpmB_rawcur_4; 
      if(name.compare("BPMA.x")==0)           val = bpmA_x; 
      if(name.compare("BPMA.y")==0)           val = bpmA_y; 
      if(name.compare("BPMB.x")==0)           val = bpmB_x; 
      if(name.compare("BPMB.y")==0)           val = bpmB_y; 
      if(name.compare("target.x")==0)         val = target_x; 
      if(name.compare("target.y")==0)         val = target_y; 
      return val;
   }
   int Print(){
      std::cout << Form("event %d, ",event) 
                << Form("run %05d, ",runNumber)
                << Form("raster 1 rawcur x = %.3lf, ",raster1_rawcur_x) 
                << Form("raster 1 rawcur y = %.3lf, ",raster1_rawcur_y) 
                << Form("raster 2 rawcur x = %.3lf, ",raster2_rawcur_x) 
                << Form("raster 2 rawcur y = %.3lf, ",raster2_rawcur_y) 
                << Form("BPMA rawcur 1 = %.3lf, ",bpmA_rawcur_1)  
                << Form("BPMA rawcur 2 = %.3lf, ",bpmA_rawcur_2)  
                << Form("BPMA rawcur 3 = %.3lf, ",bpmA_rawcur_3)  
                << Form("BPMA rawcur 4 = %.3lf, ",bpmA_rawcur_4)  
                << Form("BPMB rawcur 1 = %.3lf, ",bpmB_rawcur_1)  
                << Form("BPMB rawcur 2 = %.3lf, ",bpmB_rawcur_2)  
                << Form("BPMB rawcur 3 = %.3lf, ",bpmB_rawcur_3)  
                << Form("BPMB rawcur 4 = %.3lf, ",bpmB_rawcur_4)  
                << Form("BPMA x = %.3lf, ",bpmA_x)  
                << Form("BPMA y = %.3lf, ",bpmA_y)  
                << Form("BPMB x = %.3lf, ",bpmB_x)  
                << Form("BPMB y = %.3lf, ",bpmB_y) 
                << Form("target x = %.3lf,",target_x) 
                << Form("target y = %.3lf ",target_y) << std::endl;
      return 0; 
   } 
} beamData_t;  

// sorting rule 
bool compareBeamData_byRun(beamData s1,beamData s2){
   if(s1.runNumber < s2.runNumber) return true;
   return false;
} 

#endif 

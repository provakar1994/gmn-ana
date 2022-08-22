#ifndef UTIL_EPICS_DATA_H
#define UTIL_EPICS_DATA_H

// a simple struct for epics data 

typedef struct epicsData {
   std::string info;             // user info 
   double time;                  // EPICS time stamp        
   double IBC1H04CRCUR2;         // IBC1H04 beam current [uA] 
   double halla_p;               // beam momentum [MeV] 
   double hac_bcm_average;       // average of U1 and D1 BCM currents [uA]     
   double hac_bcm_dvm1_read;     // U1 BCM digital voltmeter readout voltage [V] 
   double hac_bcm_dvm2_read;     // D1 BCM digital voltmeter readout voltage [V]    
   double hac_bcm_dvm1_current;  // U1 BCM digital voltmeter converted to current [uA]       
   double hac_bcm_dvm2_current;  // D1 BCM digital voltmeter converted to current [uA]     
   double IPM1H04A_XPOS;         // IPM1H04A beam x position [mm] 
   double IPM1H04A_YPOS;         // IPM1H04A beam y position [mm] 
   double IPM1H04E_XPOS;         // IPM1H04E beam x position [mm] 
   double IPM1H04E_YPOS;         // IPM1H04E beam y position [mm]
   int event;                    // cumulative event number
   int runEvent;                 // run event number 
   int runNumber;                // run number 
   // constructor 
   epicsData():
      info("NONE"),time(0),IBC1H04CRCUR2(0),halla_p(0),
      hac_bcm_average(0),hac_bcm_dvm1_read(0),hac_bcm_dvm2_read(0),
      hac_bcm_dvm1_current(0),hac_bcm_dvm2_current(0),
      IPM1H04A_XPOS(0),IPM1H04A_YPOS(0),
      IPM1H04E_XPOS(0),IPM1H04E_YPOS(0),
      event(0),runEvent(0),runNumber(0)
   {}
   // get value of the member variable based on name 
   double getValue(std::string varName){
      double val=0;
      if(varName.compare("event")==0)                val = event;
      if(varName.compare("run")==0)                  val = runNumber;
      if(varName.compare("runEvent")==0)             val = runEvent;
      if(varName.compare("time")==0)                 val = time;
      if(varName.compare("halla_p")==0)              val = halla_p;
      if(varName.compare("IPM1H04A.XPOS")==0)        val = IPM1H04A_XPOS;
      if(varName.compare("IPM1H04A.YPOS")==0)        val = IPM1H04A_YPOS;
      if(varName.compare("IPM1H04E.XPOS")==0)        val = IPM1H04E_XPOS;
      if(varName.compare("IPM1H04E.YPOS")==0)        val = IPM1H04E_YPOS;
      if(varName.compare("hac_bcm_average")==0)      val = hac_bcm_average;
      if(varName.compare("hac_bcm_dvm1_read")==0)    val = hac_bcm_dvm1_read;
      if(varName.compare("hac_bcm_dvm2_read")==0)    val = hac_bcm_dvm2_read;
      if(varName.compare("hac_bcm_dvm1_current")==0) val = hac_bcm_dvm1_current;
      if(varName.compare("hac_bcm_dvm2_current")==0) val = hac_bcm_dvm2_current;
      if(varName.compare("IBC1H04CRCUR2")==0)        val = IBC1H04CRCUR2;
      return val;
   }
   // print to screen 
   int Print(){
      std::cout << Form("event %05d, "                  ,event)
	        << Form("run %05d, "                    ,runNumber)
                << Form("run event %d, "                ,runEvent)      
	        << Form("time = %.3lf, "                ,time)
	        << Form("halla_p = %.3lf, "             ,halla_p)
	        << Form("IPM1H04A_XPOS = %.3lf, "       ,IPM1H04A_XPOS)
	        << Form("IPM1H04A_YPOS = %.3lf, "       ,IPM1H04A_YPOS)
	        << Form("IPM1H04E_XPOS = %.3lf, "       ,IPM1H04E_XPOS)
	        << Form("IPM1H04E_YPOS = %.3lf, "       ,IPM1H04E_YPOS)
	        << Form("hac_bcm_average = %.3lf, "     ,hac_bcm_average)
	        << Form("hac_bcm_dvm1_read = %.3lf, "   ,hac_bcm_dvm1_read)
	        << Form("hac_bcm_dvm2_read = %.3lf, "   ,hac_bcm_dvm2_read)
	        << Form("hac_bcm_dvm1_current = %.3lf, ",hac_bcm_dvm1_current)
	        << Form("hac_bcm_dvm2_current = %.3lf, ",hac_bcm_dvm2_current)
	        << Form("IPM1H04CRCUR2 = %.3lf, "       ,IBC1H04CRCUR2) << std::endl;
      return 0;
   } 
} epicsData_t;

// sorting rule 
bool compareEPICSData_byRun(epicsData s1,epicsData s2){
   if(s1.runNumber < s2.runNumber) return true;
   return false;
} 

#endif 

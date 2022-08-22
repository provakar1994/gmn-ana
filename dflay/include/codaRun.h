#ifndef UTIL_CODA_RUN_H
#define UTIL_CODA_RUN_H

// a simple struct for CODA run info 

typedef struct codaRun {
   std::string info;       // miscellaneous user info 
   std::string rfPrefix;   // directory where ROOT files live 
   int stream;             // EVIO file stream
   int segmentBegin;       // EVIO beginning segment 
   int segmentEnd;         // EVIO end segment 
   int numFiles;           // number of ROOT files
   int runNumber;          // run number 
   // constructor 
   codaRun(): 
      info("NONE"),stream(0),segmentBegin(0),segmentEnd(0),numFiles(0),runNumber(0)
   {}
} codaRun_t; 

#endif 

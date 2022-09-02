#ifndef UTIL_CODA_RUN_H
#define UTIL_CODA_RUN_H

#include <vector>

// a simple struct for CODA run info 

typedef struct codaRun {
  std::string info;       // miscellaneous user info 
  std::string rfPrefix;   // directory where ROOT files live 
  int stream;             // EVIO file stream
  /* std::vector<int> segmentBegin;       // EVIO beginning segment  */
  /* std::vector<int> segmentEnd;         // EVIO end segment */ 
  std::vector<pair<int, int>> segmentBeg_End; // EVIO beginning & Ending segment in individual root files for a given run.
  int numFiles;           // number of ROOT files
  int runNumber;          // run number 

  // constructor 
codaRun(): 
  info("NONE"),stream(0),segmentBeg_End({}),numFiles(0),runNumber(0)
  {}

  void Clear() 
  {
    segmentBeg_End.clear();
  }

} codaRun_t; 

#endif 

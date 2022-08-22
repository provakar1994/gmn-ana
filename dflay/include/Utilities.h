#ifndef DF_UTIL_H
#define DF_UTIL_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "codaRun.h"
#include "producedVariable.h"

namespace util_df {
   unsigned long int GetUTCTimeStampFromString(std::string timeStamp,bool isDST); 
   std::string GetStringTimeStampFromUTC(unsigned long unix_time); 
   void ListDirectory(const char *path,std::vector<std::string> &list);
   int MakeDirectory(const char *path); 
   int LogMessage(const char *outpath,const char *msg,char mode='w',bool printToScreen=true); 
   int LoadRunList(const char *inpath,const char *rfPrefix,std::vector<codaRun_t> &runList); 
   int LoadBCMData(std::vector<codaRun_t> runList,BCMManager *mgr); 
   int SplitString(const char delim,const std::string myStr,std::vector<std::string> &out); 
   int GetROOTFileMetaData(const char *rfDirPath,int run,std::vector<int> &data,bool isDebug=false);
} 

#endif 

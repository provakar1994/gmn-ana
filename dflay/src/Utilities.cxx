#include "../include/Utilities.h"
//______________________________________________________________________________
namespace util_df {
   //______________________________________________________________________________
   unsigned long int GetUTCTimeStampFromString(std::string timeStamp,bool isDST){
      // takes a time stamp string and converts to UTC integer
      // returns the UTC time in SECONDS
      // note the format of the timestamp (ddd yyyy-mm-dd hh:mm::ss zzz)  
      char dayOfWeek[5],timeZone[5];
      int yy,mm,dd,hour,min,sec;
      std::sscanf(timeStamp.c_str(),"%s %04d-%02d-%02d %02d:%02d:%02d %s",
	    dayOfWeek,&yy,&mm,&dd,&hour,&min,&sec,timeZone);
      struct tm timeinfo = {0};
      timeinfo.tm_year = yy-1900;  // years since 1900 
      timeinfo.tm_mon  = mm-1;     // months since january (0--11)  
      timeinfo.tm_mday = dd;
      timeinfo.tm_hour = hour;
      timeinfo.tm_min  = min;
      timeinfo.tm_sec  = sec;
      int corr = 0;
      // FIXME: there's got to be a better way to do this... 
      // int isDST = timeinfo.tm_isdst; 
      if(isDST) corr = 3600; // daylight savings time, subtract 1 hour  
      time_t timeSinceEpoch = mktime(&timeinfo) - corr;
      return timeSinceEpoch;
   }
   //______________________________________________________________________________
   std::string GetStringTimeStampFromUTC(unsigned long unix_time){
      time_t     utime = unix_time;
      struct tm  ts;
      char       buf[100];
      // Format time as "ddd yyyy-mm-dd hh:mm:ss zzz"
      ts = *localtime(&utime);
      strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", &ts);
      std::string timeStamp = buf;
      return timeStamp;
   }
   //______________________________________________________________________________
   int MakeDirectory(const char *path){
      int rc=0;
      struct stat SB;

      const int SIZE = 200;
      char *dir_path = new char[SIZE];

      sprintf(dir_path,"%s",path);

      // check to see if the directory exists
      if( (stat(dir_path,&SB)==0)&&(S_ISDIR(SB.st_mode)) ){
	 // the directory exists!  Do nothing.  
	 // You can also check for other file types using other S_IS* macros:
	 // S_ISDIR — directory
	 // S_ISREG — regular file
	 // S_ISCHR — character device
	 // S_ISBLK — block device
	 // S_ISFIFO — FIFO
	 // S_ISLNK — symbolic link
	 // S_ISSOCK — socket
      }else{
	 rc = mkdir(dir_path,0700);
      }

      if(rc!=0) std::cout << "[MakeDirectory]: Cannot make the directory: " << path << std::endl;

      delete[] dir_path;

      return rc;
   }
   //______________________________________________________________________________
   int LoadRunList(const char *inpath,const char *rfPrefix,std::vector<codaRun_t> &runList){
      // load run list from a file 
      CSVManager *csv = new CSVManager();
      int rc = csv->ReadFile(inpath,true); // header is included 
      if(rc!=0){
	 delete csv;
	 return 1;
      }

      std::vector<std::string> data,STREAM,SEG_BEG,SEG_END,INFO;
      csv->GetColumn_byIndex_str(0,data);
      csv->GetColumn_byIndex_str(1,STREAM);
      csv->GetColumn_byIndex_str(2,SEG_BEG);
      csv->GetColumn_byIndex_str(3,SEG_END);
      csv->GetColumn_byIndex_str(4,INFO);

      codaRun_t crun;
      crun.rfPrefix = rfPrefix; // ROOT file prefix (same for all runs) 
      int aRun=0,aStr=0,aBeg=0,aEnd=0,aNumFiles=0;
      std::string info;
      const int N = data.size();
      std::vector<int> md; 
      for(int i=0;i<N;i++){
	 aRun = std::atoi( data[i].c_str() );
         rc = GetROOTFileMetaData(rfPrefix,aRun,md);
	 aStr              = md[0]; // std::atoi( STREAM[i].c_str() );
	 aBeg              = md[1]; // std::atoi( SEG_BEG[i].c_str() );
	 aEnd              = md[2]; // std::atoi( SEG_END[i].c_str() );
	 aNumFiles         = md[3]; 
	 info              = INFO[i];
	 crun.runNumber    = aRun;
	 crun.stream       = aStr;
	 crun.segmentBegin = aBeg;
	 crun.segmentEnd   = aEnd;
	 crun.numFiles     = aNumFiles;
	 crun.info         = info;
	 runList.push_back(crun);
	 // set up for next run 
	 md.clear();
      }

      delete csv;
      return 0;
   }
   //______________________________________________________________________________
   int LoadBCMData(std::vector<codaRun_t> runList,BCMManager *mgr){
      int NF=0;
      TString filePrefix,filePath;
      const int NR = runList.size();
      for(int i=0;i<NR;i++){
	 std::cout << Form("Loading run %d",runList[i].runNumber) << std::endl;
	 NF = runList[i].numFiles; 
	 if(NF==0){
	    std::cout << "No ROOT files for run " << runList[i].runNumber << std::endl;
	    continue;
	 }else{
	    std::cout << "[util_df::LoadBCMData]: Loading " << NF << " ROOT files..." << std::endl;
	    for(int j=0;j<NF;j++){
	       filePrefix = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d",
		     runList[i].rfPrefix.c_str(),runList[i].runNumber,runList[i].stream,runList[i].segmentBegin,runList[i].segmentEnd);
	       if(j==0){
		  filePath = Form("%s.root",filePrefix.Data());
	       }else{
		  filePath = Form("%s_%d.root",filePrefix.Data(),j);
	       }
	       mgr->LoadFile(filePath,runList[i].runNumber);
	    }
	    std::cout << "[util_df::LoadBCMData]: Done!" << std::endl;
	 }
      }
      return 0;
   }
   //______________________________________________________________________________
   void ListDirectory(const char *path,std::vector<std::string> &list){
      struct dirent *entry;
      DIR *dir = opendir(path);
      if(dir==NULL){
	 return;
      }
      std::string myStr;
      while ((entry = readdir(dir)) != NULL) {
	 myStr = entry->d_name;
	 list.push_back(myStr);
      }
      closedir(dir);
   }
   //______________________________________________________________________________
   int SplitString(const char delim,const std::string myStr,std::vector<std::string> &out){
      // split a string by a delimiter
      std::stringstream ss(myStr);
      std::vector<std::string> result;
      while( ss.good() ){
	 std::string substr;
	 std::getline(ss,substr,delim);
	 out.push_back(substr);
      }
      return 0;
   }
   //______________________________________________________________________________
   int LogMessage(const char *outpath,const char *msg,char mode,bool printToScreen){
      // print a log message to a file
      ios_base::openmode MODE = std::ofstream::out; 
      if(mode=='a') MODE = std::ofstream::app;  
      std::ofstream outfile;
      outfile.open(outpath,MODE);
      if( outfile.fail() ){
	 std::cout << "[util_df::LogMessage]: Cannot open the file: " << outpath << std::endl;
	 return 1;
      }else{
	 if(printToScreen) std::cout << msg << std::endl;
	 outfile << msg << std::endl;
	 outfile.close();
      } 
      return 0;  
   }
   //______________________________________________________________________________
   int GetROOTFileMetaData(const char *rfDirPath,int run,std::vector<int> &data,bool isDebug){
      // determine the beginning and end segment number for a CODA run
      // input: 
      // - rfDirPath: path to where the ROOT file(s) are located 
      // - run: CODA run number 
      // output: vector containing:  
      // - stream: stream number of the EVIO file associated with the run  
      // - begSeg: start segment number of the ROOT file associated with run 
      // - endSeg: ending segment number of the ROOT file associated with run
      // - number of files associated with the run

      int rc=-1; // assume fail 

      // first get list of files in the directory 
      std::vector<std::string> fileList; 
      ListDirectory(rfDirPath,fileList); 
      
      // skip first two entries since we get . and .. at top of list
      for(int i=0;i<2;i++) fileList.erase(fileList.begin()); 

      // identify the right index, and find number of files 
      int j=-1,fileCnt=0,theRun=0;
      const int NF = fileList.size();
      std::vector<std::string> o1;
      for(int i=0;i<NF;i++){ 
	 // split each entry based on a sufficient delimiter 
	 SplitString('_',fileList[i],o1);
	 // 3rd entry (index 2) is the one that has the run number 
	 theRun = std::atoi(o1[2].c_str());
         if(theRun==run){
	    fileCnt++; 
	    if(fileCnt==1) j = i; // found first ROOT file, save this index  
         }
	 o1.clear();
      } 

      if(fileCnt==0){
	 std::cout << "[util::GetROOTFileMetaData]: No ROOT files for run " << run << "!" << std::endl;
	 return 1;
      }

      // now parse the string indexed by j to find all the relevant info 
      std::string rfFile = fileList[j];  

      // split each entry based on a sufficient delimiter 
      std::vector<std::string> o2; 
      SplitString('_',rfFile,o1);
      // determine the stream (index 3) 
      std::string theStr = o1[3]; 
      SplitString('m',theStr,o2); 
      int stream = std::atoi(o2[1].c_str());
      o2.clear();
      // determine the segment numbers (index 4 and 5) 
      theStr = o1[4]; 
      SplitString('g',theStr,o2); 
      int begSeg = std::atoi(o2[1].c_str());  
      int endSeg = std::atoi(o1[5].c_str());

      if(isDebug){
         std::cout << Form("[Utilities::GetROOTFileMetaData]: Found run %d, stream = %d, begSeg = %d, endSeg = %d, num files = %d", 
                           run,stream,begSeg,endSeg,fileCnt) << std::endl;
      }

      data.push_back(stream); 
      data.push_back(begSeg); 
      data.push_back(endSeg); 
      data.push_back(fileCnt);
 
      return 0;
   }
} // ::util

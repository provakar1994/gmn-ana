#ifndef UTIL_ROOTFILE_MANAGER_H
#define UTIL_ROOTFILE_MANAGER_H

// a ROOT file manager

#include <cstdlib>
#include <iostream>
#include <string> 

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

// #include "Event.h"

namespace util_df {
   class ROOTFileManager { 

      private:
	 bool fIsDebug;

         std::vector< std::vector<Event<Double_t> *> > fData; // main data container 

         std::vector<std::string> fTreeName;
         std::vector<std::vector<std::string>> fBranchName; 
         std::vector<std::vector<char>> fBufSize; 

         int CheckFile(const char *filePath); 
         int LoadDataFromTree(const char *filePath,const char *treeName,std::vector<std::string> branch,std::vector<char> bufSize);  

      public: 
	 ROOTFileManager();
	 ~ROOTFileManager();

         void SetDebug(bool v=true) { fIsDebug = v; } 

         int Print(); 
	 int Clear(); 
         int PrintFileStructure(); 
	 int LoadFile(const char *filePath,const char *rfConfigPath);
	 int LoadFileStructure(const char *inpath);
         int GetTreeIndex(const char *treeName); 
         
	 TH1F *GetTH1F(const char *treeName,const char *branchName,int NBin,double min,double max);  

	 // templated methods
	 template<typename T> 
	    int GetVector(const char *treeName,const char *branchName,std::vector<T> &out){
	       // find the index corresponding to the tree name
               std::string name = treeName;
	       int k = GetTreeIndex(treeName);
	       // get the branch data
	       T val=0;
	       int NEV = fData[k].size();
               for(int i=0;i<NEV;i++){
		  val = fData[k][i]->GetData_byName(branchName);
		  out.push_back(val);
               } 
	       return 0;
	    }
 
   };  

} //::util 

#endif 

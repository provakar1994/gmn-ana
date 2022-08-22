#include "../include/ROOTFileManager.h"
//______________________________________________________________________________
namespace util_df {
   //______________________________________________________________________________
   ROOTFileManager::ROOTFileManager(){
      fIsDebug = false;
   }
   //______________________________________________________________________________
   ROOTFileManager::~ROOTFileManager(){
      int M=0;
      const int N = fData.size();
      for(int i=0;i<N;i++){
         M = fData[i].size(); 
         for(int j=0;j<M;j++) delete fData[i][j];
	 fBranchName[i].clear(); 
	 fBufSize[i].clear(); 
      } 
      fData.clear();
      fTreeName.clear();
   }
   //______________________________________________________________________________
   int ROOTFileManager::Clear(){
      // clear the data structure
      // NOTE: We are not *deleting* the vector, just setting everything 
      // back to zero effectively 
      int M=0;
      const int N = fData.size();
      for(int i=0;i<N;i++){
	 M = fData[i].size();
	 for(int j=0;j<M;j++) fData[i][j]->ClearData();
	 fBranchName[i].clear(); 
	 fBufSize[i].clear(); 
      }
      return 0;
   }
   //______________________________________________________________________________
   int ROOTFileManager::GetTreeIndex(const char *treeName){
      int k=0;
      std::string name = treeName;
      int NT = fTreeName.size();
      for(int i=0;i<NT;i++){
	 if(name.compare(fTreeName[i])==0) k = i;
      }
      return k;
   }
   //______________________________________________________________________________
   int ROOTFileManager::Print(){
      if(fIsDebug) std::cout << "[ROOTFileManager::Print] Printing data to screen..." << std::endl;
      int M=0;
      int N = fData.size();
      for(int i=0;i<N;i++){
	 std::cout << Form("Tree: %s",fTreeName[i].c_str()) << std::endl;
	 M = fData[i].size();
	 for(int j=0;j<M;j++) fData[i][j]->Print();
      }
      return 0;
   }
   //______________________________________________________________________________
   int ROOTFileManager::PrintFileStructure(){
      int NT = fTreeName.size();
      int NB = 0;
      if(NT!=0){
	 for(int i=0;i<NT;i++){
            std::cout << Form("Tree: %s",fTreeName[i].c_str()) << std::endl;
            NB = fBranchName[i].size();
            for(int j=0;j<NB;j++){
               std::cout << Form("   Branch: %s, bufSize = %c",fBranchName[i][j].c_str(),fBufSize[i][j]) << std::endl;
            }
         }
      }
      return 0;
   }
   //______________________________________________________________________________
   int ROOTFileManager::LoadFileStructure(const char *inpath){

      // check to see if we have a config loaded 
      int NB = 0;
      int NT = fTreeName.size();
      if(NT!=0){
	 std::cout << "[ROOTFileManager::LoadFileStructure]: Already have trees and branches defined!" << std::endl;
	 PrintFileStructure();
	 return 1;
      }

      // load in the list of trees and branches
      CSVManager *csv = new CSVManager();
      csv->ReadFile(inpath,true);
 
      std::vector<std::string> tree,branch,bufSize; 
      csv->GetColumn_byName_str("treeName"  ,tree   );
      csv->GetColumn_byName_str("branchName",branch );
      csv->GetColumn_byName_str("bufSize"   ,bufSize);
      delete csv; 

      // find number of branches
      std::vector<char> buf; 
      std::vector<std::string> bb;  
      std::string theTree_prev=tree[0]; 
      const int N = branch.size();
      for(int i=0;i<N;i++){
         if(tree[i].compare(theTree_prev)!=0){
            // new tree name
            fTreeName.push_back(theTree_prev);
            // store branches and buffer sizes  
            fBranchName.push_back(bb);
            fBufSize.push_back(buf);
            // clear branch history 
            bb.clear(); 
            buf.clear();
	    // need to save branch and bufSize for next tree!
	    bb.push_back(branch[i]); 
	    buf.push_back((char)bufSize[i][0]); 
         }else{
            // for a given tree name, save the branches 
            bb.push_back(branch[i]); 
            buf.push_back((char)bufSize[i][0]);  
         }        
         theTree_prev = tree[i]; 
      } 
 
      // and the last one that's not processed 
      // new tree name
      fTreeName.push_back(theTree_prev);
      // store branches and buffer sizes  
      fBranchName.push_back(bb);
      fBufSize.push_back(buf);
      // clear branch history 
      bb.clear();
      buf.clear();

      if(fIsDebug){
         std::cout << "[ROOTFileManager::LoadFileStructure]: Found " << NT << " trees: " << std::endl; 
         PrintFileStructure(); 
      }

      return 0;
   }
   //______________________________________________________________________________
   int ROOTFileManager::CheckFile(const char *filePath){
      TFile *myFile = new TFile(filePath);
      Long64_t bytesRead = myFile->GetBytesRead();
      int rc=1;
      if(bytesRead!=0){
	 // file has non-zero size
	 rc = 0;
      }
      myFile->Close();
      return rc;
   }
   //______________________________________________________________________________
   int ROOTFileManager::LoadFile(const char *filePath,const char *rfConfigPath){
      // load data from a ROOT file given the parameters defined in the rootData struct
      // data is loaded to the local CSV manager
      int rc = CheckFile(filePath);
      if(rc!=0){
	 std::cout << "[ROOTFileManager::LoadFile]: Cannot open the file " << filePath << std::endl;
	 return rc;
      }

      rc = LoadFileStructure(rfConfigPath); 
      if(rc!=0){
	 std::cout << "[ROOTFileManager::LoadFile]: Cannot open the file " << rfConfigPath << std::endl;
	 return rc;
      }

      // load data from the tree
      const int NT = fTreeName.size();
      if(fIsDebug) std::cout << "[ROOTFileManager::LoadFile]: Loading data from " << NT << " trees..." << std::endl;
      for(int i=0;i<NT;i++){
	 rc = LoadDataFromTree(filePath,fTreeName[i].c_str(),fBranchName[i],fBufSize[i]);
         if(rc!=0) std::cout << Form("[ROOTFileManager::LoadFile]: ERROR! Cannot read data for tree '%s'",fTreeName[i].c_str()) << std::endl;
      }
  
      return rc;
   }
   //______________________________________________________________________________
   int ROOTFileManager::LoadDataFromTree(const char *filePath,const char *treeName,
                                         std::vector<std::string> branch,
                                         std::vector<char> bufSize){

      if(fIsDebug) std::cout << Form("[ROOTFileManager::LoadDataFromTree]: Trying tree '%s'",treeName) << std::endl;

      // load data from a tree
      // create the TChain and attach to the tree 
      TChain *ch = nullptr;
      ch = new TChain(treeName);
      ch->Add(filePath);
      const int NN = ch->GetEntries(); // WARNING: must do this here so the tree registers in memory!

      // error checking
      if(ch==nullptr){
	 std::cout << "[ROOTFileManager::LoadDataFromTree]: ERROR! Invalid chain" << std::endl;
	 return 1;
      }

      TTree *aTree = ch->GetTree();
      if(aTree==nullptr){
	 std::cout << Form("[ROOTFileManager::LoadDataFromTree]: ERROR! Invalid tree '%s'",treeName) << std::endl;
	 return 1;
      }
       
      // passed all tests, now populate the fData container  
      const int NB = branch.size();

      if(fIsDebug) std::cout << Form("[ROOTFileManager::LoadDataFromTree]: Setting addresses for %d branches...",NB) << std::endl;

      std::vector<Int_t> var_i(NB);    // 32-bit signed integer
      std::vector<Float_t> var_f(NB);  // 32-bit floating point 
      std::vector<Double_t> var_d(NB); // 64-bit floating point 
      for(int i=0;i<NB;i++){
	 if(bufSize[i]=='D') aTree->SetBranchAddress(branch[i].c_str(),&var_d[i]);
	 if(bufSize[i]=='F') aTree->SetBranchAddress(branch[i].c_str(),&var_f[i]);
	 if(bufSize[i]=='I') aTree->SetBranchAddress(branch[i].c_str(),&var_i[i]);
      }

      if(fIsDebug) std::cout << Form("[ROOTFileManager::LoadDataFromTree]: Number of events = %d, number of branches = %d",NN,NB) << std::endl;

      std::vector< Event<Double_t> * > data; 

      for(int i=0;i<NN;i++){
	 aTree->GetEntry(i);
	 // if(fIsDebug) std::cout << Form("event %d:",i);
         // create a new event
         // FIXME: How to handle different data types? This is good for *all* branches 
         // being of type Double_t.  Be careful of casting? 
         Event<Double_t> *evt = new Event<Double_t>();
	 evt->SetID(i);  
         evt->SetVariableNames(branch);  
	 for(int j=0;j<NB;j++){ 
	    // FIXME: Technically wrong. would need unique evt classes for Double_t, Float_t, etc...
            if(bufSize[j]=='D') evt->SetData_byIndex(j,var_d[j]); 
            if(bufSize[j]=='F') evt->SetData_byIndex(j,var_f[j]); 
            if(bufSize[j]=='I') evt->SetData_byIndex(j,var_i[j]); 
         }
	 data.push_back(evt); 
      }

      fData.push_back(data); 

      if(fIsDebug) std::cout << "[ROOTFileManager::LoadDataFromTree]: Done!" << std::endl; 

      // cleanup 
      delete aTree;
      delete ch;

      return 0;
   }
   //______________________________________________________________________________
   TH1F * ROOTFileManager::GetTH1F(const char *treeName,const char *branchName,
                                   int NBin,double min,double max){
      int i   = GetTreeIndex(treeName);
      int NEV = fData[i].size();
      double arg=0;
      TString hname = Form("h%s_%s",treeName,branchName);
      TString title = Form("%s.%s" ,treeName,branchName);
      TH1F *h = new TH1F(hname,title,NBin,min,max);
      for(int j=0;j<NEV;j++){
         arg = fData[i][j]->GetData_byName(branchName); 
	 h->Fill(arg);
      }
      return h;
   }
} // ::util

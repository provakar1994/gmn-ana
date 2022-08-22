// Plot all BCM data  

#include <cstdlib>
#include <iostream>
// #include <string_view>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"
#include "TROOT.h"
#include <ROOT/RDataFrame.hxx>

#include "./include/rootData.h"
#include "./include/scalerData.h"
#include "./include/codaRun.h"
#include "./src/Event.cxx"
#include "./src/Graph.cxx"
#include "./src/CSVManager.cxx"
#include "./src/JSONManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/ROOTFileManager.cxx"

int testLogMessage(); 
int testTimeStamp(); 
int testROOTFileMetaData();
int testJSONManager(); 
int testRDataFrame();
int testROOTFileManager(); 

int Test(){

   int rc=0;
   // rc = testROOTFileMetaData();
   // rc = testJSONManager();
   // rc = testRDataFrame();
   // rc = testROOTFileManager();
   // rc = testTimeStamp();
   rc = testLogMessage();  

   return 0;
}
//______________________________________________________________________________
int testLogMessage(){
   int rc = util_df::LogMessage("test.txt","some message that I made up",'a'); 
   return rc;
}
//______________________________________________________________________________
int testTimeStamp(){

   std::string prefix   = "/lustre19/expphy/volatile/halla/sbs/flay/GMnAnalysis/rootfiles";
   std::string fileName = prefix + "/gmn_replayed-beam_13297_stream0_seg0_2.root";

   TFile *myFile = new TFile(fileName.c_str()); 

   TChain *ch = nullptr;
   ch = new TChain("T");
   ch->Add(fileName.c_str());
   int NN = ch->GetEntries();

   if(ch==nullptr){
      return 1;
   }

   TTree *aTree = ch->GetTree();
   if(aTree==nullptr) return 1;

   ULong64_t theTime=0;

   // // get the event branch 
   // TBranch *br = aTree->GetBranch("Event_Branch");
   // // get the leaf with the time stamp 
   // TLeaf *myLeaf = br->GetLeaf("fEvtHdr.fEvtTime"); 
   // myLeaf->SetAddress(&theTime); // not sure why this complains... 
   aTree->SetBranchAddress("Event_Branch.fEvtHdr.fEvtTime",&theTime);

   std::string timeStamp="";
   for(int i=0;i<100;i++){
      aTree->GetEntry(i);
      timeStamp = util_df::GetStringTimeStampFromUTC(theTime);
      std::cout << Form("event %05d, time = %llu (%s)",i,theTime,timeStamp.c_str()) << std::endl; 
   }

   return 0;
}
//______________________________________________________________________________
int testROOTFileMetaData(){

   JSONManager *jmgr = new JSONManager("./input/json/test.json");

   std::string rfPath = jmgr->GetValueFromKey_str("ROOTfile-path");
   int run=11913;
   std::vector<int> data; 
   int rc = util_df::GetROOTFileMetaData(rfPath.c_str(),run,data); 

   int stream=-1,begSeg=-1,endSeg=-1; 
   if(rc==0){
      stream = data[0]; 
      begSeg = data[1]; 
      endSeg = data[2];
      std::cout << Form("Run %d, begSeg = %d, endSeg = %d",run,begSeg,endSeg) << std::endl;
   }

   delete jmgr; 

   return 0;
}
//______________________________________________________________________________
int testJSONManager(){
   JSONManager *jmgr = new JSONManager("./input/json/test.json");
   jmgr->Print(); 

   std::string cutFile = jmgr->GetValueFromKey_str("cut-list"); 
   std::cout << cutFile << std::endl; 

   std::vector<double> cc; 
   jmgr->GetVectorFromKey<double>("calib-coeffs",cc); 

   int NCC = cc.size(); 
   for(int i=0;i<NCC;i++){
      std::cout << cc[i] << std::endl;
   }

   delete jmgr;

   return 0;
}
//______________________________________________________________________________
int testRDataFrame(){

   // path to the ROOT file
   std::string prefix   = "/lustre19/expphy/volatile/halla/sbs/flay/GMnAnalysis/rootfiles";
   std::string fileName = prefix + "/gmn_replayed-beam_13297_stream0_seg0_2.root";
 
   std::cout << "TRYING FILE: " << fileName << std::endl;

   // ROOT dataframe
   std::string_view rfName{fileName}; 
   ROOT::RDataFrame sbsDF("TSsbs",rfName,{"TSsbs"});

   // make a cut (lambda function) -- this doesn't work with tree variables with names like "branch.leaf.x"
   // auto myUnserCut = [](double sbs.bcm.unser.rate){ return sbs.bcm.unser.rate>775E+3; };   

   // make a histogram with a cut  
   // auto h = sbsDF.Filter(myUnserCut).Histo1D("sbs.bcm.unser.rate");
   auto h = sbsDF.Filter("sbs.bcm.unser.rate>775E+3").Histo1D("sbs.bcm.unser.rate");

   TCanvas *c1 = new TCanvas("c1","Test Plot",1200,800);
   c1->cd();
   h->DrawClone();
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
int testROOTFileManager(){

   std::string prefix   = "/lustre19/expphy/volatile/halla/sbs/flay/GMnAnalysis/rootfiles";
   std::string fileName = prefix + "/gmn_replayed-beam_13297_stream0_seg0_2.root";
   // path to the file that defines the ROOTfile structure
   // this is a csv of the form: treeName,branchName,bufferSize  
   // example (must include the header below)
   // treeName,branchName,bufSize 
   // TSsbs,sbs.bcm.unser.rate,D (D = Double_t, see TTree class def for details) 
   std::string structurePath = "./input/format/beam.csv";   

   util_df::ROOTFileManager *rfMgr = new util_df::ROOTFileManager();
   rfMgr->LoadFile(fileName.c_str(),structurePath.c_str());
   rfMgr->Print(); 
 
   // make a histogram  
   TH1F *h = rfMgr->GetTH1F("TSsbs","sbs.bcm.unser.rate",1000,0,900E+3);

   TCanvas *c1 = new TCanvas("c1","Test Plot",1200,800);
   c1->cd();

   h->Draw();
   c1->Update(); 

   // get a vector of data  
   std::vector<double> unserRate; 
   rfMgr->GetVector<double>("TSsbs","sbs.bcm.unser.rate",unserRate); 

   const int N = unserRate.size();
   for(int i=0;i<N;i++) std::cout << Form("%.3lf",unserRate[i]) << std::endl;
 
   delete rfMgr; 

   return 0;
}

// check the BPM data 

#include <cstdlib> 
#include <iostream>
#include <vector> 

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2D.h"

int bpmCheck(){

   gStyle->SetOptStat(0);

   int run;
   std::cout << "Enter run number: ";
   std::cin  >> run;

   TString prefix   = Form("/lustre19/expphy/volatile/halla/sbs/flay/GMnAnalysis/rootfiles");
   TString filePath = Form("%s/gmn_replayed-beam_%d_stream0_seg0_0.root",prefix.Data(),run);

   TChain *ch = new TChain("T");
   ch->Add(filePath);

   // BPMA  
   ch->Draw("Lrb.BPMA.rawcur.1:SBSrb.BPMA.rawcur.1>>ha1","","colz");
   ch->Draw("Lrb.BPMA.rawcur.2:SBSrb.BPMA.rawcur.2>>ha2","","colz");
   ch->Draw("Lrb.BPMA.rawcur.3:SBSrb.BPMA.rawcur.3>>ha3","","colz");
   ch->Draw("Lrb.BPMA.rawcur.4:SBSrb.BPMA.rawcur.4>>ha4","","colz");
   // BPMB  
   ch->Draw("Lrb.BPMB.rawcur.1:SBSrb.BPMB.rawcur.1>>hb1","","colz");
   ch->Draw("Lrb.BPMB.rawcur.2:SBSrb.BPMB.rawcur.2>>hb2","","colz");
   ch->Draw("Lrb.BPMB.rawcur.3:SBSrb.BPMB.rawcur.3>>hb3","","colz");
   ch->Draw("Lrb.BPMB.rawcur.4:SBSrb.BPMB.rawcur.4>>hb4","","colz");

   const int N = 4;
   TH2D **ha = new TH2D*[N];
   TH2D **hb = new TH2D*[N];
   for(int i=0;i<N;i++){
      ha[i] = (TH2D *)gDirectory->Get( Form("ha%d",i+1) );
      hb[i] = (TH2D *)gDirectory->Get( Form("hb%d",i+1) );
   }

   TCanvas *c2 = new TCanvas("c2",Form("Run %d: BPMA",run),0,0,1400,1000);
   c2->Divide(2,2);

   for(int i=0;i<N;i++){
      c2->cd(i+1);
      ha[i]->Draw();
      c2->Update();
   }

   TCanvas *c3 = new TCanvas("c3",Form("Run %d: BPMB",run),0,0,1400,1000);
   c3->Divide(2,2);

   for(int i=0;i<N;i++){
      c3->cd(i+1);
      hb[i]->Draw();
      c3->Update();
   }

   return 0;
} 

// check the raster data 

#include <cstdlib> 
#include <iostream>
#include <vector> 

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2D.h"

int rasterCheck(){

   gStyle->SetOptStat(0);

   int run;
   std::cout << "Enter run number: ";
   std::cin  >> run;

   TString prefix   = Form("/lustre19/expphy/volatile/halla/sbs/flay/GMnAnalysis/rootfiles");
   TString filePath = Form("%s/gmn_replayed-beam_%d_stream0_seg0_0.root",prefix.Data(),run);

   TChain *ch = new TChain("T");
   ch->Add(filePath);

   ch->Draw("Lrb.Raster.rawcur.x:SBSrb.Raster.rawcur.x>>hxx","","colz");
   ch->Draw("Lrb.Raster.rawcur.y:SBSrb.Raster.rawcur.y>>hyy","","colz");
   ch->Draw("Lrb.Raster2.rawcur.x:SBSrb.Raster2.rawcur.x>>h2xx","","colz");
   ch->Draw("Lrb.Raster2.rawcur.y:SBSrb.Raster2.rawcur.y>>h2yy","","colz");

   TH2D *hxx  = (TH2D *)gDirectory->Get("hxx");
   TH2D *hyy  = (TH2D *)gDirectory->Get("hyy");
   TH2D *h2xx = (TH2D *)gDirectory->Get("h2xx");
   TH2D *h2yy = (TH2D *)gDirectory->Get("h2yy");

   TCanvas *c2 = new TCanvas("c2","Raster Size Analysis",0,0,1400,1000);
   c2->Divide(2,2);

   c2->cd(1);
   hxx->GetXaxis()->SetTitle("SBSrb.Raster.rawcur.x");
   hxx->GetXaxis()->CenterTitle();
   hxx->GetYaxis()->SetTitle("Lrb.Raster.rawcur.x");
   hxx->GetYaxis()->CenterTitle();
   hxx->Draw("colz");
   c2->Update();

   c2->cd(2);
   hyy->GetXaxis()->SetTitle("SBSrb.Raster.rawcur.y");
   hyy->GetXaxis()->CenterTitle();
   hyy->GetYaxis()->SetTitle("Lrb.Raster.rawcur.y");
   hyy->GetYaxis()->CenterTitle();
   hyy->Draw("colz");
   c2->Update();

   c2->cd(3);
   h2xx->GetXaxis()->SetTitle("SBSrb.Raster2.rawcur.x");
   h2xx->GetXaxis()->CenterTitle();
   h2xx->GetYaxis()->SetTitle("Lrb.Raster2.rawcur.x");
   h2xx->GetYaxis()->CenterTitle();
   h2xx->Draw("colz");
   c2->Update();

   c2->cd(4);
   h2yy->GetXaxis()->SetTitle("SBSrb.Raster2.rawcur.y");
   h2yy->GetXaxis()->CenterTitle();
   h2yy->GetYaxis()->SetTitle("Lrb.Raster2.rawcur.y");
   h2yy->GetYaxis()->CenterTitle();
   h2yy->Draw("colz");
   c2->Update();

   return 0;
} 

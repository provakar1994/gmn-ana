// read in analysis results and make comparison plots 

#include <cstdlib>
#include <iostream>

#include "./include/bcm.h"
#include "./include/codaRun.h"
#include "./include/charge.h"
#include "./src/ABA.cxx"
#include "./src/CSVManager.cxx"
#include "./src/BCMManager.cxx"
#include "./src/Utilities.cxx"
#include "./src/JSONManager.cxx"
#include "./src/Graph.cxx"
#include "./src/bcmUtilities.cxx"
#include "./src/cutUtilities.cxx"

int GetCC(std::string bcmName,std::vector<calibCoeff_t> data,std::vector<double> &x,std::vector<double> &dx); 

int bcmCompareResults(const char *confPath){
   
   int rc=0;

   JSONManager *jmgr = new JSONManager(confPath);
   std::vector<std::string> confList;
   jmgr->GetVectorFromKey_str("configs",confList); 
   delete jmgr;  

   const int N = 7; 
   std::string bcm[N] = {"u1","unew","d1","d3","d10","dnew"};

   std::vector<calibCoeff_t> cc; 

   int M=0;
   char inpath[200]; 

   std::vector<double> x,dx;
   std::vector<double> conf;
   std::vector<double> off_uns,offErr_uns,ped_uns,pedErr_uns,slope_uns,slopeErr_uns; 
   std::vector<double> off_u1,offErr_u1,ped_u1,pedErr_u1,slope_u1,slopeErr_u1; 
   std::vector<double> off_unew,offErr_unew,ped_unew,pedErr_unew,slope_unew,slopeErr_unew; 
   std::vector<double> off_d1,offErr_d1,ped_d1,pedErr_d1,slope_d1,slopeErr_d1; 
   std::vector<double> off_d3,offErr_d3,ped_d3,pedErr_d3,slope_d3,slopeErr_d3; 
   std::vector<double> off_d10,offErr_d10,ped_d10,pedErr_d10,slope_d10,slopeErr_d10; 
   std::vector<double> off_dnew,offErr_dnew,ped_dnew,pedErr_dnew,slope_dnew,slopeErr_dnew; 

   // loop over configs
   const int NCONF = confList.size(); 
   for(int i=0;i<NCONF;i++){
      conf.push_back( (double)(i+1) );
      std::cout << "Loading " << confList[i] << std::endl;
      sprintf(inpath,"./output/%s/result.csv",confList[i].c_str()); 
      rc = bcm_util::LoadCalibrationCoefficients(inpath,cc);
      // get data for each BCM
      M = cc.size();
      for(int j=0;j<M;j++){
	 GetCC(bcm[j],cc,x,dx);
         if(bcm[j].compare("unser")==0){
	    ped_uns.push_back(x[0]); 
	    off_uns.push_back(x[1]); 
	    slope_uns.push_back(x[2]); 
	    pedErr_uns.push_back(dx[0]); 
	    offErr_uns.push_back(dx[1]); 
	    slopeErr_uns.push_back(dx[2]); 
         }else if(bcm[j].compare("u1")==0){
	    ped_u1.push_back(x[0]); 
	    off_u1.push_back(x[1]); 
	    slope_u1.push_back(x[2]); 
	    pedErr_u1.push_back(dx[0]); 
	    offErr_u1.push_back(dx[1]); 
	    slopeErr_u1.push_back(dx[2]); 
         }else if(bcm[j].compare("unew")==0){
	    ped_unew.push_back(x[0]); 
	    off_unew.push_back(x[1]); 
	    slope_unew.push_back(x[2]); 
	    pedErr_unew.push_back(dx[0]); 
	    offErr_unew.push_back(dx[1]); 
	    slopeErr_unew.push_back(dx[2]); 
         }else if(bcm[j].compare("d1")==0){
	    ped_d1.push_back(x[0]); 
	    off_d1.push_back(x[1]); 
	    slope_d1.push_back(x[2]); 
	    pedErr_d1.push_back(dx[0]); 
	    offErr_d1.push_back(dx[1]); 
	    slopeErr_d1.push_back(dx[2]); 
         }else if(bcm[j].compare("d3")==0){
	    ped_d3.push_back(x[0]); 
	    off_d3.push_back(x[1]); 
	    slope_d3.push_back(x[2]); 
	    pedErr_d3.push_back(dx[0]); 
	    offErr_d3.push_back(dx[1]); 
	    slopeErr_d3.push_back(dx[2]); 
         }else if(bcm[j].compare("d10")==0){
	    ped_d10.push_back(x[0]); 
	    off_d10.push_back(x[1]); 
	    slope_d10.push_back(x[2]); 
	    pedErr_d10.push_back(dx[0]); 
	    offErr_d10.push_back(dx[1]); 
	    slopeErr_d10.push_back(dx[2]); 
         }else if(bcm[j].compare("dnew")==0){
	    ped_dnew.push_back(x[0]); 
	    off_dnew.push_back(x[1]); 
	    slope_dnew.push_back(x[2]); 
	    pedErr_dnew.push_back(dx[0]); 
	    offErr_dnew.push_back(dx[1]); 
	    slopeErr_dnew.push_back(dx[2]); 
         }
	 // clear for next bcm 
	 x.clear();
	 dx.clear();
      }
      // clear for next config 
      cc.clear();
   }
 
   // create plots
   // TGraphErrors *gOff_uns  = graph_df::GetTGraphErrors(conf,off_uns ,offErr_uns);  
   TGraphErrors *gOff_u1   = graph_df::GetTGraphErrors(conf,off_u1  ,offErr_u1);  
   TGraphErrors *gOff_unew = graph_df::GetTGraphErrors(conf,off_unew,offErr_unew);  
   TGraphErrors *gOff_d1   = graph_df::GetTGraphErrors(conf,off_d1  ,offErr_d1);  
   TGraphErrors *gOff_d3   = graph_df::GetTGraphErrors(conf,off_d3  ,offErr_d3);  
   TGraphErrors *gOff_d10  = graph_df::GetTGraphErrors(conf,off_d10 ,offErr_d10);  
   TGraphErrors *gOff_dnew = graph_df::GetTGraphErrors(conf,off_dnew,offErr_dnew);  

   TGraphErrors *gPed_u1   = graph_df::GetTGraphErrors(conf,ped_u1  ,pedErr_u1);  
   TGraphErrors *gPed_unew = graph_df::GetTGraphErrors(conf,ped_unew,pedErr_unew);  
   TGraphErrors *gPed_d1   = graph_df::GetTGraphErrors(conf,ped_d1  ,pedErr_d1);  
   TGraphErrors *gPed_d3   = graph_df::GetTGraphErrors(conf,ped_d3  ,pedErr_d3);  
   TGraphErrors *gPed_d10  = graph_df::GetTGraphErrors(conf,ped_d10 ,pedErr_d10);  
   TGraphErrors *gPed_dnew = graph_df::GetTGraphErrors(conf,ped_dnew,pedErr_dnew);  

   TGraphErrors *gSlope_u1   = graph_df::GetTGraphErrors(conf,slope_u1  ,slopeErr_u1);  
   TGraphErrors *gSlope_unew = graph_df::GetTGraphErrors(conf,slope_unew,slopeErr_unew);  
   TGraphErrors *gSlope_d1   = graph_df::GetTGraphErrors(conf,slope_d1  ,slopeErr_d1);  
   TGraphErrors *gSlope_d3   = graph_df::GetTGraphErrors(conf,slope_d3  ,slopeErr_d3);  
   TGraphErrors *gSlope_d10  = graph_df::GetTGraphErrors(conf,slope_d10 ,slopeErr_d10);  
   TGraphErrors *gSlope_dnew = graph_df::GetTGraphErrors(conf,slope_dnew,slopeErr_dnew);  

   graph_df::SetParameters(gPed_u1,20,kBlack); 
   graph_df::SetParameters(gOff_u1,20,kBlack); 
   graph_df::SetParameters(gSlope_u1,20,kBlack); 

   graph_df::SetParameters(gPed_unew,20,kBlack); 
   graph_df::SetParameters(gOff_unew,20,kBlack); 
   graph_df::SetParameters(gSlope_unew,20,kBlack); 

   graph_df::SetParameters(gPed_d1,20,kBlack); 
   graph_df::SetParameters(gOff_d1,20,kBlack); 
   graph_df::SetParameters(gSlope_d1,20,kBlack); 

   graph_df::SetParameters(gPed_d3,20,kBlack); 
   graph_df::SetParameters(gOff_d3,20,kBlack); 
   graph_df::SetParameters(gSlope_d3,20,kBlack);
 
   graph_df::SetParameters(gPed_d10,20,kBlack); 
   graph_df::SetParameters(gOff_d10,20,kBlack); 
   graph_df::SetParameters(gSlope_d10,20,kBlack); 

   graph_df::SetParameters(gPed_dnew,20,kBlack); 
   graph_df::SetParameters(gOff_dnew,20,kBlack); 
   graph_df::SetParameters(gSlope_dnew,20,kBlack); 

   double xSize = 0.05;
   double ySize = 0.06; 

   TCanvas *c1 = new TCanvas("c1","Upstream BCM",1200,600);
   c1->Divide(3,2);

   c1->cd(1);
   gPed_u1->Draw("ap");
   graph_df::SetLabels(gPed_u1,"u1 Pedestal [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gPed_u1,xSize,ySize); 
   gPed_u1->Draw("ap");
   c1->Update(); 

   c1->cd(2);
   gOff_u1->Draw("ap");
   graph_df::SetLabels(gOff_u1,"u1 Offset [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gOff_u1,xSize,ySize); 
   gOff_u1->Draw("ap");
   c1->Update(); 

   c1->cd(3);
   gSlope_u1->Draw("ap");
   graph_df::SetLabels(gSlope_u1,"u1 Gain [Hz/#muA]","Analysis",""); 
   graph_df::SetLabelSizes(gSlope_u1,xSize,ySize); 
   gSlope_u1->Draw("ap");
   c1->Update(); 

   c1->cd(4);
   gPed_unew->Draw("ap");
   graph_df::SetLabels(gPed_unew,"unew Pedestal [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gPed_unew,xSize,ySize); 
   gPed_unew->Draw("ap");
   c1->Update(); 

   c1->cd(5);
   gOff_unew->Draw("ap");
   graph_df::SetLabels(gOff_unew,"unew Offset [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gOff_unew,xSize,ySize); 
   gOff_unew->Draw("ap");
   c1->Update(); 

   c1->cd(6);
   gSlope_unew->Draw("ap");
   graph_df::SetLabels(gSlope_unew,"unew Gain [Hz/#muA]","Analysis",""); 
   graph_df::SetLabelSizes(gSlope_unew,xSize,ySize); 
   gSlope_unew->Draw("ap");
   c1->Update(); 

   TCanvas *c2 = new TCanvas("c2","Downstream BCM: d1 and d3",1200,600);
   c2->Divide(3,2);

   c2->cd(1);
   gPed_d1->Draw("ap");
   graph_df::SetLabels(gPed_d1,"d1 Pedestal [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gPed_d1,xSize,ySize); 
   gPed_d1->Draw("ap");
   c2->Update(); 

   c2->cd(2);
   gOff_d1->Draw("ap");
   graph_df::SetLabels(gOff_d1,"d1 Offset [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gOff_d1,xSize,ySize); 
   gOff_d1->Draw("ap");
   c2->Update(); 

   c2->cd(3);
   gSlope_d1->Draw("ap");
   graph_df::SetLabels(gSlope_d1,"d1 Gain [Hz/#muA]","Analysis",""); 
   graph_df::SetLabelSizes(gSlope_d1,xSize,ySize); 
   gSlope_d1->Draw("ap");
   c2->Update(); 

   c2->cd(4);
   gPed_d3->Draw("ap");
   graph_df::SetLabels(gPed_d3,"d3 Pedestal [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gPed_d3,xSize,ySize); 
   gPed_d3->Draw("ap");
   c2->Update(); 

   c2->cd(5);
   gOff_d3->Draw("ap");
   graph_df::SetLabels(gOff_d3,"d3 Offset [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gOff_d3,xSize,ySize); 
   gOff_d3->Draw("ap");
   c2->Update(); 

   c2->cd(6);
   gSlope_d3->Draw("ap");
   graph_df::SetLabels(gSlope_d3,"d3 Gain [Hz/#muA]","Analysis",""); 
   graph_df::SetLabelSizes(gSlope_d3,xSize,ySize); 
   gSlope_d3->Draw("ap");
   c2->Update(); 

   TCanvas *c3 = new TCanvas("c3","Downstream BCM: d10 and dnew",1200,600);
   c3->Divide(3,2);

   c3->cd(1);
   gPed_d10->Draw("ap");
   graph_df::SetLabels(gPed_d10,"d10 Pedestal [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gPed_d10,xSize,ySize); 
   gPed_d10->Draw("ap");
   c3->Update(); 

   c3->cd(2);
   gOff_d10->Draw("ap");
   graph_df::SetLabels(gOff_d10,"d10 Offset [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gOff_d10,xSize,ySize); 
   gOff_d10->Draw("ap");
   c3->Update(); 

   c3->cd(3);
   gSlope_d10->Draw("ap");
   graph_df::SetLabels(gSlope_d10,"d10 Gain [Hz/#muA]","Analysis",""); 
   graph_df::SetLabelSizes(gSlope_d10,xSize,ySize); 
   gSlope_d10->Draw("ap");
   c3->Update(); 

   c3->cd(4);
   gPed_dnew->Draw("ap");
   graph_df::SetLabels(gPed_dnew,"dnew Pedestal [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gPed_dnew,xSize,ySize); 
   gPed_dnew->Draw("ap");
   c3->Update(); 

   c3->cd(5);
   gOff_dnew->Draw("ap");
   graph_df::SetLabels(gOff_dnew,"dnew Offset [Hz]","Analysis",""); 
   graph_df::SetLabelSizes(gOff_dnew,xSize,ySize); 
   gOff_dnew->Draw("ap");
   c3->Update(); 

   c3->cd(6);
   gSlope_dnew->Draw("ap");
   graph_df::SetLabels(gSlope_dnew,"dnew Gain [Hz/#muA]","Analysis",""); 
   graph_df::SetLabelSizes(gSlope_dnew,xSize,ySize); 
   gSlope_dnew->Draw("ap");
   c3->Update(); 

   return 0;
}
//_____________________________________________________________________________
int GetCC(std::string bcmName,std::vector<calibCoeff_t> data,std::vector<double> &x,std::vector<double> &dx){
   // get calibration coefficient values for bcmName
   const int N = data.size();
   for(int i=0;i<N;i++){
      if(data[i].dev.compare(bcmName)==0){
	 x.push_back(data[i].pedestal); 
	 x.push_back(data[i].offset); 
	 x.push_back(data[i].slope); 
	 dx.push_back(data[i].pedestalErr); 
	 dx.push_back(data[i].offsetErr); 
	 dx.push_back(data[i].slopeErr); 
      }
   }
   return 0;
}

/*
  This script will compare the new vs old pedestal subtraction algorithms.
  Old algo: just looked at the first 4 FADC time bins to get the average 
  pedestal. New algo: Looks at first as well as last 4 bins and takes the
  minimum between the two to get average pedestal ("4" is user-configurable).
  -------
  P. Datta <pdbforce@jlab.org> Created July 26, 2022
*/

#include <iostream>
#include <vector>
#include <TH1D.h>
#include <TFile.h>
#include <TEventList.h>

#include "../../include/Constants.h"
#include "../../src/SetROOTVar.cpp"
#include "../../src/ExpConstants.cpp"

void pedAlgo_study(const char* rfile="",TCut cut="",
		   const char* outFile="pedAlgo_study_sbs11_nalgo.root")
{
  
  TChain *C = new TChain("T");
  C->Add(rfile);
  cout << " Opened ROOT files with " << C->GetEntries() << " events" << endl;
  TEventList *elist = new TEventList("elist");
  C->Draw(">>elist",cut);
  cout << " No of events passed the cut = " << elist->GetN() << endl;
  
  int max=1000;
  C->SetBranchStatus("*",0);

  //bbsh variables
  int aSHndata;
  double aSH[max],apSH[max],acSH[max],ampSH[max],amppSH[max],elSH[max],eSH,idblkSH;
  std::vector<std::string> shvar = {"a","a","a_p","a_c","a_amp","a_amp_p","adcelemID","e","idblk"};
  std::vector<void*> shvar_mem = {&aSH,&aSHndata,&apSH,&acSH,&ampSH,&amppSH,&elSH,&eSH,&idblkSH};
  setrootvar::setbranch(C,"bb.sh",shvar,shvar_mem,1);

  //bbps variables
  int aPSndata;
  double aPS[max],apPS[max],acPS[max],ampPS[max],amppPS[max],elPS[max],ePS,idblkPS;
  std::vector<std::string> psvar = {"a","a","a_p","a_c","a_amp","a_amp_p","adcelemID","e","idblk"};
  std::vector<void*> psvar_mem = {&aPS,&aPSndata,&apPS,&acPS,&ampPS,&amppPS,&elPS,&ePS,&idblkPS};
  setrootvar::setbranch(C,"bb.ps",psvar,psvar_mem,1);

  //hcal variables
  int aHCALndata;
  double aHCAL[max],apHCAL[max],acHCAL[max],ampHCAL[max],amppHCAL[max],
    elHCAL[max],eHCAL,idblkHCAL;
  std::vector<std::string> hcalvar = {"a","a","a_p","a_c","a_amp","a_amp_p",
				      "adcelemID","e","idblk"};
  std::vector<void*> hcalvar_mem = {&aHCAL,&aHCALndata,&apHCAL,&acHCAL,&ampHCAL,
				  &amppHCAL,&elHCAL,&eHCAL,&idblkHCAL};
  setrootvar::setbranch(C,"sbs.hcal",hcalvar,hcalvar_mem,1);

  //MC variables
  double mc_sigma, mc_omega, mc_np, mc_ep, mc_fnucl;
  std::vector<std::string> mc = {"mc_sigma","mc_omega","mc_np","mc_ep","mc_fnucl"};
  std::vector<void*> mc_mem = {&mc_sigma,&mc_omega,&mc_np,&mc_ep,&mc_fnucl};
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  //TFile *fout = new TFile(outFile,"RECREATE");
  TString data_file = "hcal_ML_short.csv";
  ofstream ml_data;
  ml_data.open(data_file);
  ml_data << "Event No., True nucleon momentum (GeV/c), n(0) or p(1), ADC, (pedestal subtracted), (pC), for, all, 288, channels, ...," << endl;

  TH1D *h_e_SH = new TH1D("h_e_SH","CLuster energy (GeV)",300,0,5);
  TH1D *h_a_SH = new TH1D("h_a_SH","ADC integral(pC) with Pedestal",600,0,600);
  TH1D *h_a_p_SH = new TH1D("h_a_p_SH","ADC integral(pC)",400,0,400);
  TH1D *h_a_c_SH = new TH1D("h_a_c_SH","ADC integral(GeV)",200,0,4);
  TH1D *h_a_amp_SH = new TH1D("h_a_amp_SH","ADC amplitude(mV) with Pedestal",800,0,800);
  TH1D *h_a_amp_p_SH = new TH1D("h_a_amp_p_SH","ADC amplitude(mV)",600,0,600);
  
  TH1D *h_e_PS = new TH1D("h_e_PS","CLuster energy (GeV)",300,0,4);
  TH1D *h_a_PS = new TH1D("h_a_PS","ADC integral(pC) with Pedestal",300,0,300);
  TH1D *h_a_p_PS = new TH1D("h_a_p_PS","ADC integral(pC)",300,0,300);
  TH1D *h_a_c_PS = new TH1D("h_a_c_PS","ADC integral(GeV)",200,0,4);
  TH1D *h_a_amp_PS = new TH1D("h_a_amp_PS","ADC amplitude(mV) with Pedestal",600,0,600);
  TH1D *h_a_amp_p_PS = new TH1D("h_a_amp_p_PS","ADC amplitude(mV)",400,0,400);

  double temp[288];
  for(int i=0; i<288; i++){
    temp[i] = 0.0;      
  }

  long nevent=0,nevents=elist->GetN();
  while( C->GetEntry(elist->GetEntry(nevent++)) ){
    if( nevent % 1000 == 0 ) cout << nevent << "/" << nevents << "\r";
    cout.flush();
    
    // Shower 
    for(int ihit=0;ihit<aSHndata;ihit++){
      // if((int)elSH[ihit]==105){
      // if((int)idblkSH==105){
	h_e_SH->Fill( eSH );
	h_a_SH->Fill( aSH[ihit] );
	h_a_p_SH->Fill( apSH[ihit] );
	h_a_c_SH->Fill( acSH[ihit] );
	h_a_amp_SH->Fill( ampSH[ihit] );
	h_a_amp_p_SH->Fill( amppSH[ihit] );   
	//}
    }

    // PreShwoer
    for(int ihit=0;ihit<aPSndata;ihit++){
      // if((int)elPS[ihit]==25){
      // if((int)idblkPS==25){
	h_e_PS->Fill( ePS );
	h_a_PS->Fill( aPS[ihit] );
	h_a_p_PS->Fill( apPS[ihit] );
	h_a_c_PS->Fill( acPS[ihit] );
	h_a_amp_PS->Fill( ampPS[ihit] );
	//h_a_amp_p_PS->Fill( amppPS[ihit] );
	//}
    }

    ml_data << nevent << "," << mc_fnucl << "," << mc_np << ",";
    // HCAL
    for(int ihit=0;ihit<aHCALndata;ihit++){
      temp[int(elHCAL[ihit])] = apHCAL[ihit];
    }

    for(int i=0; i<288; i++){
      ml_data << temp[i] << ",";   
      temp[i] = 0.0;
    }
    
    ml_data << endl;
  }

  ml_data.close();

  //fout->Write();
  elist->Delete();
}

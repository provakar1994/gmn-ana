/* 
   This macro will determine the n and p detection efficiencies of HCAL using g4sbs's particle
   gun generator. We digitize and replay the simulated data to make it more realistic.
   -----
   Config files: sbs4-sbs0p-simu/conf_hcal_det_effi_p(n)gun.json
   -----
   P. Datta  Created  09-12-2022 (Based on AJR Puckett's momentum calibration script)
*/
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../../include/gmn-ana.h"
#include "../../../dflay/src/JSONManager.cxx"

double pi = constant::pi;
double Mp = constant::Mp;
double Mn = constant::Mn;

TH1D* MakeHisto(int, int, double, double, std::string);

int hcal_det_effi_gun(const char *configfilename, 
		      std::string filebase="siout/hcal_det_effi")
{

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfiles = jmgr->GetValueFromKey_str("rootfiles");
  TChain *C = new TChain("T");
  C->Add(rootfiles.c_str());

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // setting up global cuts
  // std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  // TCut globalcut = gcut.c_str();
  // TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // hcal clus var
  double eHCAL, xHCAL, yHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  //MC variables
  double mc_vx, mc_vy, mc_vz, mc_ep, mc_epx, mc_epy, mc_epz, mc_sigma, mc_omega;
  std::vector<std::string> mc = {"mc_vx","mc_vy","mc_vz","mc_ep","mc_epx","mc_epy","mc_epz","mc_sigma","mc_omega"};
  std::vector<void*> mc_mem = {&mc_vx,&mc_vy,&mc_vz,&mc_ep,&mc_epx,&mc_epy,&mc_epz,&mc_sigma,&mc_omega};
  setrootvar::setbranch(C, "MC", mc, mc_mem);

  // defining the outputfile
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");
  std::string outFile = Form("%s_%s_gun_sbs%d_sbs%dp.root", filebase.c_str(), Ntype.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag());
  TFile *fout = new TFile(outFile.c_str(), "RECREATE");

  // defining histograms
  TH1D *h_dxHCAL = new TH1D("h_dxHCAL",Form("%s gun; x_{HCAL} - x_{exp} (m);",Ntype.c_str()), 250, -2.5, 2.5);
  TH1D *h_dyHCAL = new TH1D("h_dyHCAL",Form("%s gun; y_{HCAL} - y_{exp} (m);",Ntype.c_str()), 250, -1.25, 1.25);
  TH2F *h2_xy = new TH2F("h2_xy",Form("%s gun",Ntype.c_str()),12,-0.92837,0.92837,24,-2.35183,1.45182);
  TH2F *h2_xy_exp = new TH2F("h2_xy_exp",Form("%s gun",Ntype.c_str()),12,-0.92837,0.92837,24,-2.35183,1.45182);

  std::vector<double> h_eHCAL_lims; jmgr->GetVectorFromKey<double>("h_eHCAL_lims",h_eHCAL_lims);
  TH1D *h_eHCAL = new TH1D("h_eHCAL",Form("%s gun; HCAL Cluster Energy (GeV);",Ntype.c_str()),int(h_eHCAL_lims[0]),h_eHCAL_lims[1],h_eHCAL_lims[2]);

  std::vector<double> h_pN_lims; jmgr->GetVectorFromKey<double>("h_pN_lims",h_pN_lims);
  TH2F *h2_eHCAL_pN = new TH2F("h2_eHCAL_pN",Form("%s gun;Nucleon Momentum (GeV/c);HCAL Cluster Energy (GeV)",Ntype.c_str()),int(h_pN_lims[0]),h_pN_lims[1],h_pN_lims[2],int(h_eHCAL_lims[0]),h_eHCAL_lims[1],h_eHCAL_lims[2]);
  TProfile *h2_eHCAL_pN_prof = new TProfile("h2_eHCAL_pN_prof",Form("%s gun;Nucleon Momentum (GeV/c);HCAL Cluster Energy (GeV)",Ntype.c_str()),int(h_pN_lims[0]),h_pN_lims[1],h_pN_lims[2],h_eHCAL_lims[1],h_eHCAL_lims[2]);
  TH2F *h2_eHCAL_pN_cut = new TH2F("h2_eHCAL_pN_cut",Form("%s gun | Threshold = E_{mean} / 4. per bin;Nucleon Momentum (GeV/c);HCAL Cluster Energy (GeV)",Ntype.c_str()),int(h_pN_lims[0]),h_pN_lims[1],h_pN_lims[2],int(h_eHCAL_lims[0]),h_eHCAL_lims[1],h_eHCAL_lims[2]);

  TH1D *h_effi = new TH1D("h_effi",Form("%s gun; Nucleon Momentum (GeV/c); Efficiency (p)",Ntype.c_str()),int(h_pN_lims[0]),h_pN_lims[1],h_pN_lims[2]);

  TTree *Tout = new TTree("Tout", "Active area cut implemented before filling");
  double T_pN;         Tout->Branch("pN", &T_pN, "pN/D");                      // nucleon momentum
  double T_eHCAL;      Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D");             // HCAL cluster energy  
  double T_xHCAL;      Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D");             // observed x (from HCAL cluster)
  double T_yHCAL;      Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D");             // observed y (from HCAL cluster)
  double T_xHCAL_exp;  Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); // expected x (from MC truth)
  double T_yHCAL_exp;  Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D"); // expected y (from MC truth)

  // costruct axes of HCAL CoS in Hall CoS
  vector<TVector3> HCAL_axes;
  kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  double HCALheight_offset = jmgr->GetValueFromKey<double>("HCALheight_offset");
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + HCALheight_offset*HCAL_axes[0];  

  // define the histograms
  TH1D *h_eHCAL_array[int(h_pN_lims[0])];
  TH1D *h_eHCAL_array_cut[int(h_pN_lims[0])];
  for (int ibin=0; ibin<int(h_pN_lims[0]); ibin++) {
    double hmax = h_eHCAL_lims[2];
    // variable ranges for splitted histograms
    if (ibin < int(h_pN_lims[0])/1.6) hmax = h_eHCAL_lims[2] / 1.2;
    if (ibin < int(h_pN_lims[0])/4) hmax = h_eHCAL_lims[2] / 2.5;
    if (ibin < int(h_pN_lims[0])/8) hmax = h_eHCAL_lims[2] / 4.0;
    h_eHCAL_array[ibin] = MakeHisto(ibin, h_eHCAL_lims[0], h_eHCAL_lims[1], hmax, "");
    h_eHCAL_array_cut[ibin] = MakeHisto(ibin, h_eHCAL_lims[0], h_eHCAL_lims[1], hmax, "_cut");
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 500); c1->cd();

  // looping through the tree ---------------------------------------
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0;
  while (C->GetEntry(nevent++)) {
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // apply global cuts efficiently (AJRP method)
    // currenttreenum = C->GetTreeNumber();
    // if (nevent == 1 || currenttreenum != treenum) {
    //   treenum = currenttreenum;
    //   GlobalCut->UpdateFormulaLeaves();
    // } 
    // bool passedgCut = GlobalCut->EvalInstance(0) != 0;
    // if (passedgCut) {

    // constructing q-vector
    TVector3 vertex(mc_vx, mc_vy, mc_vz);
    double thetaN_expect = acos(mc_epz / mc_ep); // final state N's theta
    double phiN_expect = atan2(mc_epy, mc_epx);  // final state N's phi
    TVector3 pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);

    // HCAL active area cut
    vector<double> hcal_active_area = cut::hcal_active_area_simu(); // Exc. 1 blk from all 4 sides
    if (!cut::inHCAL_activeA(xyHCAL_exp[0], xyHCAL_exp[1], hcal_active_area)) continue; 
    // if (!cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area)) continue; // Introduces bias!!   

    T_pN = mc_ep;
    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;
    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    Tout->Fill();

    h_dxHCAL->Fill(xHCAL - xyHCAL_exp[0]);
    h_dyHCAL->Fill(yHCAL - xyHCAL_exp[1]);
    h2_xy->Fill(yHCAL, xHCAL);
    h2_xy_exp->Fill(xyHCAL_exp[1], xyHCAL_exp[0]); 

    h_eHCAL->Fill(eHCAL);
    h2_eHCAL_pN->Fill(mc_ep, eHCAL);   
    h2_eHCAL_pN_prof->Fill(mc_ep, eHCAL, 1.);

    // filling histo for diff. pN values
    double binwidth = (h_pN_lims[2] - h_pN_lims[1]) / h_pN_lims[0];
    int ibin = (mc_ep - h_pN_lims[1]) / binwidth; // h2_eHCAL_pN->GetFixBin(eHCAL); should work too!
    h_eHCAL_array[ibin]->Fill(eHCAL);

    //} // globalcuts   
  } // while event loop
  std::cout << std::endl;

  h2_eHCAL_pN_prof->SetMarkerStyle(7);
  h2_eHCAL_pN_prof->SetMarkerColor(2);

  // Determine bin wise threshold (E_peak / 4.)
  cout << endl << "Determining thresholds.." << endl;
  vector<double> threshold;
  TF1 *fgaus = new TF1("fgaus", "gaus");
  for (int ibin=0; ibin<int(h_pN_lims[0]); ibin++) {
    // let's fit the histograms to get E_mean (or, E_peak)
    int maxBin = h_eHCAL_array[ibin]->GetMaximumBin();
    double maxBinCenter = h_eHCAL_array[ibin]->GetXaxis()->GetBinCenter( maxBin );
    double maxCount = h_eHCAL_array[ibin]->GetMaximum();
    double binWidth = h_eHCAL_array[ibin]->GetBinWidth(maxBin);
    double stdDev = h_eHCAL_array[ibin]->GetStdDev();

    // Reject low energy peak
    if (maxBin < 2) {
      while ( h_eHCAL_array[ibin]->GetBinContent(maxBin+1) < h_eHCAL_array[ibin]->GetBinContent(maxBin) || 
	      h_eHCAL_array[ibin]->GetBinContent(maxBin+1) == h_eHCAL_array[ibin]->GetBinContent(maxBin) ) {
	maxBin++;
      }
      h_eHCAL_array[ibin]->GetXaxis()->SetRange(maxBin+1, h_eHCAL_array[ibin]->GetNbinsX());
      maxBin = h_eHCAL_array[ibin]->GetMaximumBin();
      maxBinCenter = h_eHCAL_array[ibin]->GetXaxis()->GetBinCenter(maxBin);
      maxCount = h_eHCAL_array[ibin]->GetMaximum();
      binWidth = h_eHCAL_array[ibin]->GetBinWidth(maxBin);
      stdDev = h_eHCAL_array[ibin]->GetStdDev();
    }

    // let's fit
    double lowerBinC = h_eHCAL_lims[1] + maxBin*binWidth - stdDev;
    double upperBinC = h_eHCAL_lims[1] + maxBin*binWidth + stdDev;
    fgaus->SetParameters(maxCount, maxBinCenter, stdDev);
    fgaus->SetRange(lowerBinC, upperBinC);
    h_eHCAL_array[ibin]->Fit(fgaus, "RQ");
    threshold.push_back(fgaus->GetParameter(1) / 4.);  // threshold = E_peak/4
  }

  // 2nd event loop
  cout << endl << "Looping over data again to apply threshold cuts.." << endl;
  nevent = 0;
  while (C->GetEntry(nevent++)) {
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // constructing q-vector
    TVector3 vertex(mc_vx, mc_vy, mc_vz);
    double thetaN_expect = acos(mc_epz / mc_ep); // final state N's theta
    double phiN_expect = atan2(mc_epy, mc_epx);  // final state N's phi
    TVector3 pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);

    // HCAL active area cut
    vector<double> hcal_active_area = cut::hcal_active_area_simu(); // Exc. 1 blk from all 4 sides
    if (!cut::inHCAL_activeA(xyHCAL_exp[0], xyHCAL_exp[1], hcal_active_area)) continue; 

    // filling histo for diff. pN values
    double binwidth = (h_pN_lims[2] - h_pN_lims[1]) / h_pN_lims[0];
    int ibin = (mc_ep - h_pN_lims[1]) / binwidth; // h2_eHCAL_pN->GetFixBin(eHCAL); should work too!
    if (eHCAL > threshold[ibin]) {  // threshold cut
      h_eHCAL_array_cut[ibin]->Fill(eHCAL);
      double binwidthy = (h_eHCAL_lims[2] - h_eHCAL_lims[1]) / h_eHCAL_lims[0];
      int ibiny = (eHCAL - h_eHCAL_lims[1]) / binwidthy;
      h2_eHCAL_pN_cut->SetBinContent(ibin, ibiny, 1);
    }
  }
  
  cout << endl << endl << "Calculating detection efficiencies.." << endl;
  for (int ibin=0; ibin<int(h_pN_lims[0]); ibin++) {
    double efficiency = (h_eHCAL_array_cut[ibin]->Integral() / h_eHCAL_array[ibin]->Integral()) * 100.;
    double binwidth = (h_pN_lims[2] - h_pN_lims[1]) / h_pN_lims[0];
    double binx = h_pN_lims[1] + ibin*binwidth + binwidth/2.;
    h_effi->Fill(binx, efficiency);
  }

  h_effi->SetLineWidth(0);
  h_effi->SetMarkerStyle(7);
  if (Ntype.compare("p")) h_effi->SetMarkerColor(2);
  else h_effi->SetMarkerColor(1);
  h_effi->GetYaxis()->SetRangeUser(80.,105.);
  h_effi->Draw();

  fout->Write();

  cout << endl << "------" << endl;
  cout << " Outfile : " << outFile << endl;
  cout << "------" << endl;

  delete jmgr;
  return 0;
}


// ---------------- Create generic histogram function ----------------
TH1D* MakeHisto(int ibin, int bins, double min, double max, std::string suf="")
{
  TH1D *h = new TH1D(TString::Format("h_eHCAL%s_%d", suf.c_str(), ibin),
		     TString::Format("h_eHCAL%s_%d", suf.c_str(), ibin), bins, min, max);
  // h->SetStats(0);
  // h->SetLineWidth(2);
  // h->GetYaxis()->SetLabelSize(0.1);
  // h->GetYaxis()->SetLabelOffset(-0.17);
  // h->GetYaxis()->SetNdivisions(5);
  return h;
}

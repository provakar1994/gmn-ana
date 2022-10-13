/* 
   This macro will determine the n and p detection efficiencies of HCAL using g4sbs's particle
   gun generator. We digitize and replay the simulated data to make it more realistic.
   -----
   P. Datta  Created  09-12-2022 (Based on AJR Puckett's momentum calibration script)
*/
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../../include/gmn-ana.h"
#include "../../../dflay/src/JSONManager.cxx"

double pi = constant::pi;
double Mp = constant::Mp;
double Mn = constant::Mn;

int hcal_det_effi_pgun(const char *configfilename, 
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

  // parameters for normalization
  std::string targetType = jmgr->GetValueFromKey_str("target_type");
  double I_beam = jmgr->GetValueFromKey<double>("beam_current");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

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
  TH1F *h_eHCAL = new TH1F("h_eHCAL","; HCAL Cluster Energy (GeV);", 200, 0., 2.);
  TH1F *h_eHCAL_cut = new TH1F("h_eHCAL_cut","; HCAL Cluster Energy w Cut (GeV);", 200, 0., 2.);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);", 250, -2.5, 2.5);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25, 1.25);
  TH2F *h_xy = new TH2F("h_xy", "",12,-0.92837,0.92837,24,-2.35183,1.45182);
  TH2F *h_xy_cut = new TH2F("h_xy_cut", "",12,-0.92837,0.92837,24,-2.35183,1.45182);
  TH2F *h_xy_exp = new TH2F("h_xy_exp", "",12,-0.92837,0.92837,24,-2.35183,1.45182);
  TH2F *h_xy_exp_cut = new TH2F("h_xy_exp_cut", "",12,-0.92837,0.92837,24,-2.35183,1.45182);

  TTree *Tout = new TTree("Tout", "Active area cut implemented before filling");
  double T_eHCAL;      Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D");
  double T_xHCAL;      Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D");
  double T_yHCAL;      Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D");
  double T_xHCAL_exp;  Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D");
  double T_yHCAL_exp;  Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");

  // costruct axes of HCAL CoS in Hall CoS
  vector<TVector3> HCAL_axes;
  kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  double HCALheight_offset = jmgr->GetValueFromKey<double>("HCALheight_offset");
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + HCALheight_offset*HCAL_axes[0];  

  // looping through the tree ---------------------------------------
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0;
  while (C->GetEntry(nevent++)) {
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (passedgCut) {

      // constructing q-vector
      TVector3 vertex(mc_vx, mc_vy, mc_vz);
      double thetaN_expect = acos(mc_epz / mc_ep); // final state N's theta
      double phiN_expect = atan2(mc_epy, mc_epx);  // final state N's phi
      TVector3 pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);

      // Expected position of the q vector at HCAL
      vector<double> xyHCAL_exp; 
      kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);

      // HCAL active area cut
      vector<double> hcal_active_area = cut::hcal_active_area_simu(); // Exc. 1 blk from all 4 sides
      // if (!cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area)) continue; // Introduces bias!!   
      if (!cut::inHCAL_activeA(xyHCAL_exp[0], xyHCAL_exp[1], hcal_active_area)) continue; 

      T_eHCAL = eHCAL;
      T_xHCAL = xHCAL;
      T_yHCAL = yHCAL;
      T_xHCAL_exp = xyHCAL_exp[0];
      T_yHCAL_exp = xyHCAL_exp[1];
      Tout->Fill();

      h_dxHCAL->Fill(xHCAL - xyHCAL_exp[0]);
      h_dyHCAL->Fill(yHCAL - xyHCAL_exp[1]);
      h_xy->Fill(yHCAL, xHCAL);
      h_xy_exp->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);    

      h_eHCAL->Fill(eHCAL);
      if (eHCAL < 0.03775) {
	h_eHCAL_cut->Fill(eHCAL);
	h_xy_cut->Fill(yHCAL, xHCAL);
	h_xy_exp_cut->Fill(xyHCAL_exp[1], xyHCAL_exp[0]); 
      }

    } // globalcuts   
  } // while event loop
  std::cout << std::endl;

  fout->Write();

  cout << endl << "------" << endl;
  cout << " Outfile : " << outFile << endl;
  cout << "------" << endl;

  delete jmgr;
  return 0;
}

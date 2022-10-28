/* 
   This macro will determine nucleon detection efficiency of HCAL using beam data.
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

#include "../../include/gmn-ana.h"
#include "../../dflay/src/JSONManager.cxx"

double pi = constant::pi;
double Mp = constant::Mp;
double Mn = constant::Mn;

int hcal_det_effi_data (const char *configfilename, std::string filebase="pdout/hcal_det_effi")
{

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  //std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  TChain *C = new TChain("T");
  C->Add(rootfile_dir.c_str());
  // for (int i=0; i<runnums.size(); i++) {
  //   std::string rfname = rootfile_dir + Form("/*%d*",runnums[i]);
  //   C->Add(rfname.c_str());
  // }

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  int model = jmgr->GetValueFromKey<int>("model");
  if ((model != 0) && (model != 1)) {
    std::cerr << "Enter a valid model number! **!**" << std::endl;
    throw;
  }

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // hcal clus var
  double eHCAL, xHCAL, yHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // tdctrig variable (N/A for simulation)
  int tdcElemN;
  double tdcTrig[maxNtr], tdcElem[maxNtr];
  std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.ps.e", 1);
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);

  // defining the outputfile
  std::string outFile = Form("%s_sbs%d_sbs%dp_model%d_data.root", filebase.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), model);
  TFile *fout = new TFile(outFile.c_str(), "RECREATE");

  // defining histograms
  TH1D *h_W = new TH1D("h_W",";W (GeV);", 250, 0., 2.);
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1F *h_Q2 = new TH1F("h_Q2", "Q^{2} distribution", 100, 0., 5.);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);", 250, -2.5, 2.5);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25, 1.25);

  // Do the energy loss calculation here ...........

  // costruct axes of HCAL CoS in Hall CoS
  vector<TVector3> HCAL_axes;
  kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + kine::HCALOriginOffset(HCAL_axes, "data");  

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
    if (!passedgCut) continue;

    // coin time cut (N/A for simulation)
    double bbcal_time=0., hcal_time=0.;
    for(int ihit=0; ihit<tdcElemN; ihit++){
      if(tdcElem[ihit]==5) bbcal_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_time=tdcTrig[ihit];
    }
    double coin_time = hcal_time - bbcal_time;
    //h_coin_time->Fill( coin_time );
    if(fabs(coin_time - 510) > 10.) continue;

    // kinematic parameters
    double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV) [Get it from EPICS, eventually]
    double ebeam_corr = ebeam; //- MeanEloss;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    /* Reaction    : e + e' -> N + N' */
    /* Conservation: Pe + Peprime = PN + PNprime */
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);                 
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN("p", PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, "p");

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + pi; 
    /* Different modes of calculation. Goal is to achieve the best resolution
       model 0 = uses reconstructed p as independent variable
       model 1 = uses reconstructed angles as independent variable */
    if (model == 0) {
      nu = ebeam_corr - Peprime.E();
      pN_expect = kine::pN_expect(nu, "p");
      thetaN_expect = acos((ebeam_corr - Peprime.Pz()) / pN_expect);
    } else if (model == 1) {
      nu = ebeam_corr - pcentral;
      pN_expect = kine::pN_expect(nu, "p");
      thetaN_expect = acos((ebeam_corr - pcentral*cos(etheta)) / pN_expect);
    }
    TVector3 pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);


    double Q2recon = kine::Q2(ebeam_corr, Peprime.E(), etheta);
    h_Q2->Fill(Q2recon); 
    double Wrecon = kine::W(ebeam_corr, Peprime.E(), Q2recon, "p");
    h_W->Fill(Wrecon);
    double dpel = Peprime.E()/pcentral - 1.0;
    h_dpel->Fill(dpel);

    // W cut
    if (abs(Wrecon - 0.88) > 0.4) continue;

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    h_dxHCAL->Fill(xHCAL - xyHCAL_exp[0]);
    h_dyHCAL->Fill(yHCAL - xyHCAL_exp[1]);

    //HCAL active area and safety margin cuts
    vector<double> hcal_active_area = cut::hcal_active_area_data(); // Exc. 1 blk from all 4 sides
    // vector<double> hcal_safety_margin = cut::hcal_safety_margin(delx_sigma, dely_sigma, hcal_active_area);
    // if (inHCAL_activeA(xHCAL, yHCAL, hcal_active_area)); //do stuff
    // if (inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], delx_shift, hcal_safety_margin)); //do stuff

      
  } // while event loop
  std::cout << std::endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Divide(2,1);
  c1->cd(1); h_dxHCAL->Draw();
  c1->cd(2); h_dyHCAL->Draw();
  // c1->cd(1); h_W->Draw();
  // c1->cd(2); h_dpel->Draw();

  fout->Write();

  delete jmgr;
  return 0;
}

/* 
   This macro will be used to get elastic yields for gmn.
   E.g. Config. File: sbs14-sbs70p/conf_qelas_ana_data.json
   -----
   P. Datta  Created  11-02-2022 
*/
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gmn-ana.h"
#include "../../dflay/src/JSONManager.cxx"

int qelas_ana_data (const char *configfilename, std::string filebase="pdout/qelas_ana_data")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  TChain *C = new TChain("T");
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    std::string rfname = rootfile_dir + Form("/*%d*",runnums[i]);
    C->Add(rfname.c_str());
  }

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = jmgr->GetValueFromKey<int>("model");
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

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
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL};
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
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d_data.root", 
			 filebase.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), model);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  //temp
  TH1D *h_dpel_p = new TH1D("h_dpel_p",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  TH1D *h_dpel_n = new TH1D("h_dpel_n",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim; jmgr->GetVectorFromKey<double>("h_dxHCAL_lims", hdx_lim);
  vector<double> hdy_lim; jmgr->GetVectorFromKey<double>("h_dyHCAL_lims", hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_fcut = new TH1F("h_dxHCAL_fcut","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_afcut = new TH1F("h_dxHCAL_afcut","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  // TH1F *h_coin_time = new TH1F("h_coin_time", "Coincidence time (ns)", 200, 380, 660);

  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // Do the energy loss calculation here ...........

  // HCAL cut definitions
  double sbs_kick = jmgr->GetValueFromKey<double>("sbs_kick");
  vector<double> dx_p; jmgr->GetVectorFromKey<double>("dx_p", dx_p);
  vector<double> dy_p; jmgr->GetVectorFromKey<double>("dy_p", dy_p);
  double Nsigma_cut_dx_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_p");
  double Nsigma_cut_dy_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_p");
  vector<double> dx_n; jmgr->GetVectorFromKey<double>("dx_n", dx_n);
  vector<double> dy_n; jmgr->GetVectorFromKey<double>("dy_n", dy_n);
  double Nsigma_cut_dx_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_n");
  double Nsigma_cut_dy_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_n");
  vector<double> hcal_active_area = cut::hcal_active_area_data(); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // elastic cut limits
  double Wmin = jmgr->GetValueFromKey<double>("Wmin");
  double Wmax = jmgr->GetValueFromKey<double>("Wmax");

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromKey<double>("hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromKey<double>("hcal_hoffset");
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
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

    // coin time cut (N/A for simulation)  !! Not a reliable cut - loosing a lot of elastics
    // double bbcal_time=0., hcal_time=0.;
    // for(int ihit=0; ihit<tdcElemN; ihit++){
    //   if(tdcElem[ihit]==5) bbcal_time=tdcTrig[ihit];
    //   if(tdcElem[ihit]==0) hcal_time=tdcTrig[ihit];
    // }
    // double coin_time = hcal_time - bbcal_time;  h_coin_time->Fill( coin_time );
    // if(fabs(coin_time - 518) > 10.) continue;

    // kinematic parameters
    double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV) [Get it from EPICS, eventually]
    double ebeam_corr = ebeam; //- MeanEloss;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    /* Reaction    : e + e' -> N + N'
       Conservation: Pe + Peprime = PN + PNprime */
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);                 
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN(Ntype, PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype);

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi; 
    /* Different modes of calculation. Goal is to achieve the best resolution
       model 0 = uses reconstructed p as independent variable
       model 1 = uses reconstructed angles as independent variable 
       model 2 => uses 4-vector calculation */
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon, W2recon;
    if (model == 0) {
      nu = ebeam_corr - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((ebeam_corr - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(ebeam_corr, Peprime.E(), etheta);
      W2recon = kine::W2(ebeam_corr, Peprime.E(), Q2recon, Ntype);
    } else if (model == 1) {
      nu = ebeam_corr - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((ebeam_corr - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(ebeam_corr, Peprime.E(), etheta);
      W2recon = kine::W2(ebeam_corr, Peprime.E(), Q2recon, Ntype);
    } else if (model == 2) {
      TLorentzVector q = Pe - Peprime; // 4-momentum of virtual photon
      pNhat = q.Vect().Unit();
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
    }

    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0; h_dpel->Fill(dpel);

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL - xyHCAL_exp[0];  
    double dy = yHCAL - xyHCAL_exp[1]; 

    if (Wrecon > Wmin && Wrecon < Wmax && abs(dy) < 0.3) h_dxHCAL->Fill(dx);

    // HCAL active area and safety margin cuts [Fiducial region]
    bool AR_cut = cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);
    // if (!cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area)) continue;
    // if (!cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_magin)) continue: 

    // W cut
    //if (Wrecon < Wmin || Wrecon > Wmax) continue;
    //if (abs(Wrecon - 0.876) > 0.2) continue;
    h_W->Fill(Wrecon);
    if (Wrecon > Wmin && Wrecon < Wmax) {

      // histos with fiducial cut
      if (AR_cut && FR_cut) {
	if (abs(dy) < 0.3) h_dxHCAL_fcut->Fill(dx);
	h_dyHCAL->Fill(dy);
	h2_rcHCAL->Fill(cblkHCAL, rblkHCAL);
	h2_dxdyHCAL->Fill(dy, dx);
      } else {
	if (abs(dy) < 0.3) h_dxHCAL_afcut->Fill(dx);
      }
    }

    // HCAL cut
    bool pcut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;
    bool ncut = pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) + pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2) <= 1.;
    if (pcut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
    if (ncut) h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
    if (pcut || ncut) { 
      h_W_cut->Fill(Wrecon);
    } else {
      h_W_acut->Fill(Wrecon);
    }
      
  } // event loop
  std::cout << std::endl << std::endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1); h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);
 
  c1->cd(2); // h2_rcHCAL->Draw("colz");
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3); //h_dxHCAL->Draw();
  h2_xyHCAL_p->Draw("colz");
  TLine L1h_AR;
  L1h_AR.SetLineColor(2); L1h_AR.SetLineWidth(4); L1h_AR.SetLineStyle(9);
  L1h_AR.DrawLine(hcal_active_area[2],hcal_active_area[1],hcal_active_area[3],hcal_active_area[1]);
  TLine L2h_AR;
  L2h_AR.SetLineColor(2); L2h_AR.SetLineWidth(4); L2h_AR.SetLineStyle(9);
  L2h_AR.DrawLine(hcal_active_area[2],hcal_active_area[0],hcal_active_area[3],hcal_active_area[0]);
  TLine L1v_AR;
  L1v_AR.SetLineColor(2); L1v_AR.SetLineWidth(4); L1v_AR.SetLineStyle(9);
  L1v_AR.DrawLine(hcal_active_area[2],hcal_active_area[0],hcal_active_area[2],hcal_active_area[1]);
  TLine L2v_AR;
  L2v_AR.SetLineColor(2); L2v_AR.SetLineWidth(4); L2v_AR.SetLineStyle(9);
  L2v_AR.DrawLine(hcal_active_area[3],hcal_active_area[0],hcal_active_area[3],hcal_active_area[1]);
  // safety margin
  TLine L1h_SM;
  L1h_SM.SetLineColor(4); L1h_SM.SetLineWidth(4); L1h_SM.SetLineStyle(9);
  L1h_SM.DrawLine(hcal_safety_margin[2],hcal_safety_margin[1],hcal_safety_margin[3],hcal_safety_margin[1]);
  TLine L2h_SM;
  L2h_SM.SetLineColor(4); L2h_SM.SetLineWidth(4); L2h_SM.SetLineStyle(9);
  L2h_SM.DrawLine(hcal_safety_margin[2],hcal_safety_margin[0],hcal_safety_margin[3],hcal_safety_margin[0]);
  TLine L1v_SM;
  L1v_SM.SetLineColor(4); L1v_SM.SetLineWidth(4); L1v_SM.SetLineStyle(9);
  L1v_SM.DrawLine(hcal_safety_margin[2],hcal_safety_margin[0],hcal_safety_margin[2],hcal_safety_margin[1]);
  TLine L2v_SM;
  L2v_SM.SetLineColor(4); L2v_SM.SetLineWidth(4); L2v_SM.SetLineStyle(9);
  L2v_SM.DrawLine(hcal_safety_margin[3],hcal_safety_margin[0],hcal_safety_margin[3],hcal_safety_margin[1]);

  c1->cd(4); 
  h2_xyHCAL_n->Draw("colz");
  //TLine L1h_AR;
  L1h_AR.SetLineColor(2); L1h_AR.SetLineWidth(4); L1h_AR.SetLineStyle(9);
  L1h_AR.DrawLine(hcal_active_area[2],hcal_active_area[1],hcal_active_area[3],hcal_active_area[1]);
  //TLine L2h_AR;
  L2h_AR.SetLineColor(2); L2h_AR.SetLineWidth(4); L2h_AR.SetLineStyle(9);
  L2h_AR.DrawLine(hcal_active_area[2],hcal_active_area[0],hcal_active_area[3],hcal_active_area[0]);
  //TLine L1v_AR;
  L1v_AR.SetLineColor(2); L1v_AR.SetLineWidth(4); L1v_AR.SetLineStyle(9);
  L1v_AR.DrawLine(hcal_active_area[2],hcal_active_area[0],hcal_active_area[2],hcal_active_area[1]);
  //TLine L2v_AR;
  L2v_AR.SetLineColor(2); L2v_AR.SetLineWidth(4); L2v_AR.SetLineStyle(9);
  L2v_AR.DrawLine(hcal_active_area[3],hcal_active_area[0],hcal_active_area[3],hcal_active_area[1]);
  // safety margin
  //TLine L1h_SM;
  L1h_SM.SetLineColor(4); L1h_SM.SetLineWidth(4); L1h_SM.SetLineStyle(9);
  L1h_SM.DrawLine(hcal_safety_margin[2],hcal_safety_margin[1],hcal_safety_margin[3],hcal_safety_margin[1]);
  //TLine L2h_SM;
  L2h_SM.SetLineColor(4); L2h_SM.SetLineWidth(4); L2h_SM.SetLineStyle(9);
  L2h_SM.DrawLine(hcal_safety_margin[2],hcal_safety_margin[0],hcal_safety_margin[3],hcal_safety_margin[0]);
  //TLine L1v_SM;
  L1v_SM.SetLineColor(4); L1v_SM.SetLineWidth(4); L1v_SM.SetLineStyle(9);
  L1v_SM.DrawLine(hcal_safety_margin[2],hcal_safety_margin[0],hcal_safety_margin[2],hcal_safety_margin[1]);
  //TLine L2v_SM;
  L2v_SM.SetLineColor(4); L2v_SM.SetLineWidth(4); L2v_SM.SetLineStyle(9);
  L2v_SM.DrawLine(hcal_safety_margin[3],hcal_safety_margin[0],hcal_safety_margin[3],hcal_safety_margin[1]);

  outFile.ReplaceAll(".root",".png");
  c1->Print(outFile.Data(),"png");

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  fout->Write();
  sw->Delete();
  delete jmgr;
  return 0;
}

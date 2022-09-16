/* 
   This macro will determine the n and p detection efficiencies of HCAL using 
   simulated data. We digitize and replay the simulated data to make it more realistic.
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

#include "../../../include/Constants.h"  // namespace constant
#include "../../../src/SetROOTVar.cpp"   // namespace setrootvar
#include "../../../src/ExpConstants.cpp" // namespace expconst
#include "../../../src/KinematicVar.cpp" // namespace kine

#include "../../../dflay/src/JSONManager.cxx"

double pi = constant::pi;
double Mp = constant::Mp;
double Mn = constant::Mn;

int hcal_det_effi_simu(const char *configfilename, 
		       std::string filebase="hcal_det_effi")
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
  SBSconfig sbsconf;
  expconst::LoadSBSconfig(conf, sbsmag, sbsconf);
  sbsconf.Print();

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  int model = jmgr->GetValueFromKey<int>("model");
  if ((model != 0) && (model != 1)) {
    std::cerr << "Enter a valid model number! **!**" << std::endl;
    throw;
  }

  // parameters for normalization
  std::string targetType = jmgr->GetValueFromKey_str("target_type");
  double I_beam = jmgr->GetValueFromKey<double>("beam_current");
  // total generated simulation events per job
  double ngen_total = jmgr->GetValueFromKey<double>("ngen_total"); 
  double lumi = kine::Luminosity(I_beam, targetType);

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy (N/A for simulation)
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // bbsh clus var
  double eSH, xSH, ySH;
  std::vector<std::string> shclvar = {"e","x","y"};
  std::vector<void*> shclvar_mem = {&eSH,&xSH,&ySH};
  setrootvar::setbranch(C, "bb.sh", shclvar, shclvar_mem);

  // bbps clus var
  double ePS, xPS, yPS;
  std::vector<std::string> psclvar = {"e","x","y"};
  std::vector<void*> psclvar_mem = {&ePS,&xPS,&yPS};
  setrootvar::setbranch(C, "bb.ps", psclvar, psclvar_mem);

  // hcal clus var
  double eHCAL, xHCAL, yHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz",
  				    "r_x","r_y","r_th","r_ph",
  				    "tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,
  				  &xfp,&yfp,&thfp,&phfp,
  				  &xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C, "bb.tr", trvar, trvar_mem);

  // tdctrig variable (N/A for simulation)
  // int tdcElemN;
  // double tdcTrig[maxNtr], tdcElem[maxNtr];
  // std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  // std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  // setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);

  //MC variables
  double mc_sigma, mc_omega;
  std::vector<std::string> mc = {"mc_sigma","mc_omega"};
  std::vector<void*> mc_mem = {&mc_sigma,&mc_omega};
  setrootvar::setbranch(C, "MC", mc, mc_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);

  // defining the outputfile
  std::string outFile = Form("%s_sbs%d_sbs%dp_model%d_simu.root", 
			     filebase.c_str(), sbsconf.sbsconf, sbsconf.sbsmag, model);
  TFile *fout = new TFile(outFile.c_str(), "RECREATE");

  // defining histograms
  TH1D *h_W = new TH1D("h_W",";W (GeV);", 250, 0., 2.);
  TH1F *h_Q2 = new TH1F("h_Q2", "Q^{2} distribution", 100, 0., 5.);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);", 250, -2.5, 2.5);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25, 1.25);

  // Do the energy loss calculation here ...........

  // costruct axes of HCAL CoS in Hall CoS
  vector<TVector3> HCAL_axes;
  kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  double HCALheight_offset = jmgr->GetValueFromKey<double>("HCALheight_offset");
  TVector3 HCAL_origin = sbsconf.HCALdist*HCAL_axes[2] + HCALheight_offset*HCAL_axes[0];  

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
      // cross section weighted normalization factor
      double weight = mc_sigma*mc_omega*lumi / ngen_total;

      // kinematic parameters
      double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV)
      /* HALLA_p N/A for simulation */
      // double ebeam_epics = HALLA_p / 1E-3;     // Beam energy from Tree (GeV)
      // if (ebeam_epics > 1) {                  // use value from tree when reasonable
      // 	if (abs(ebeam_epics - ebeam) > 0.1) {
      // 	  std::cout 
      // 		    << " *!* WARNING: > 100 MeV fluctuation in beam energy!" << std::endl
      // 		    << " event " << nevent << " ebeam_expect = " << ebeam 
      // 		    << " HALLA_p = " << ebeam_epics << std::endl << " ----- " << std::endl;
      // 	}
      // 	ebeam = ebeam_epics;
      // }
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
      double Ntheta_expect = 0.;        // expected recoil nucleon theta
      double Nphi_expect = ephi + pi; 
      /* Different modes of calculation. Goal is to achieve the best resolution
	 model 0 = uses reconstructed p as independent variable
	 model 1 = uses reconstructed angles as independent variable */
      if (model == 0) {
	nu = ebeam_corr - Peprime.E();
	pN_expect = kine::pN_expect(nu, "p");
	Ntheta_expect = acos((ebeam_corr - Peprime.Pz()) / pN_expect);
      } else if (model == 1) {
	nu = ebeam_corr - pcentral;
	pN_expect = kine::pN_expect(nu, "p");
	Ntheta_expect = acos((ebeam_corr - pcentral*cos(etheta)) / pN_expect);
      }
      TVector3 pNhat = kine::qVect_unit(Ntheta_expect, Nphi_expect);

      // Expected position of the q vector at HCAL
      vector<double> xyHCAL_exp; 
      kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
      h_dxHCAL->Fill(xHCAL - xyHCAL_exp[0]);
      h_dyHCAL->Fill(yHCAL - xyHCAL_exp[1]);

      double Q2recon = kine::Q2(ebeam_corr, Peprime.E(), etheta);
      h_Q2->Fill(Q2recon); 
      double Wrecon = kine::W(ebeam_corr, Peprime.E(), Q2recon, "p");
      h_W->Fill(Wrecon);
      double dpel = Peprime.E()/pcentral - 1.0;

    } // if passedgCut
    
  } // while event loop
  std::cout << std::endl;

  fout->Write();
  
  delete jmgr;
  return 0;
}

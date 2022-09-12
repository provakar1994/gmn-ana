/* 
   This macro will determine the n and p detection efficiencies of HCAL using 
   simulated data. We digitize and replay the simulated data to make it more realistic.
   -----
   P. Datta  Created  09-12-2022 (Based on AJR Puckett's momentum calibration script)
*/
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TMath.h"
#include "TChain.h"
#include "TLatex.h"
#include "TVector3.h"
#include "TRotation.h"

#include "../../../include/Constants.h"
#include "../../../src/SetROOTVar.cpp"
#include "../../../src/ExpConstants.cpp"

#include "../../../dflay/src/JSONManager.cxx"

double pi = constant.pi;
double Mp = constant.Mp;
double Mn = constant.Mn;

int hcal_detection_effi(const char *configfilename, const char *outputfilename="hcal_dxdy_ld2.root")
{

  // reading input config file
  JSONManager *jmgr = new JSONManager(confPath);
  TCut globalcut = jmgr->GetValueFromKey_str("global_cut");
  
  // seting up the desired SBS configuration
  SBSconfig sbsconf;
  expconst::LoadSBSconfig(conf, sbsconf);
  sbsconf.Print();

  
  
  delete jmgr;

  return 0;
}

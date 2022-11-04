#ifndef PD_UTIL_H
#define PD_UTIL_H

#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"

#include "../include/ExpConstants.h"

namespace util_pd {

  /* #################################################
     ##                HCAL Histograms              ##  
     ################################################# */
  TH2F *TH2FHCALface_rc(std::string name);      // returns TH2F for HCAL face (row,col)
  TH2F *TH2FHCALface_xy_data(std::string name); // returns TH2F for HCAL face (x,y) [Data]
  TH2F *TH2FHCALface_xy_simu(std::string name); // returns TH2F for HCAL face (x,y) [Simu]
  TH2F *TH2FdxdyHCAL(std::string name);   // returns TH2F for dxdyHCAL

  /* #################################################
     ##              Kinematic Histograms           ##  
     ################################################# */
  TH1F *TH1FhW(std::string name);   // returns W histogram
  TH1F *TH1FhQ2(std::string name,   // returns Q2 histogram
		int conf);   // SBS config
}

#endif

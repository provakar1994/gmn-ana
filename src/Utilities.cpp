#include "../include/Utilities.h"

namespace util_pd {

  /* #################################################
     ##                HCAL Related                 ##  
     ################################################# */
  //_____________________________________
  TH2F *TH2FHCALface_rc(std::string name) {
    // returns TH2F for HCAL face (row,col)
    /* NOTE: HCAL block id (ibblk) starts from 1 and goes up to 288 but both
             HCAL row (rowblk) and column starts from 0 and goes up to 23 and 
             11, respectively. Extremely annoying! */
    TH2F *h = new TH2F(name.c_str(), ";HCAL columns;HCAL rows",
		       expconst::hcalcol, 0, expconst::hcalcol,
		       expconst::hcalrow, 0, expconst::hcalrow);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_data(std::string name) {
    // returns TH2F for HCAL face (x,y) [Data]
    double y_min = expconst::yHCAL_r_DB - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB + expconst::hcalblk_h/2.;
    TH2F *h = new TH2F(name.c_str(), ";y_{HCAL} (m);x_{HCAL} (m)",
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_simu(std::string name) {
    // returns TH2F for HCAL face (x,y) [Simu]
    double y_min = expconst::yHCAL_r_DB_MC - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB_MC + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB_MC - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB_MC + expconst::hcalblk_h/2.;
    TH2F *h = new TH2F(name.c_str(), ";y_{HCAL} (m);x_{HCAL} (m)",
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  TH2F *TH2FdxdyHCAL(std::string name) {
    // returns TH2F for dxdyHCAL
    TH2F *h = new TH2F(name.c_str(), "; y_{HCAL} - y_{exp} (m); x_{HCAL} - x_{exp} (m)",
		       250, -1.25, 1.25, 250, -3.5, 2);
    return h;
  }

  /* #################################################
     ##              Kinematic Histograms           ##  
     ################################################# */
  TH1F *TH1FhW(std::string name) {
    // returns W histogram
    TH1F *h = new TH1F(name.c_str(), "W Distribution (GeV)", 250,0,2);
    return h;
  }
  TH1F *TH1FhQ2(std::string name,       // Name of histogram
		int conf) {             // SBS config
    // returns Q2 histogram
    int nbin; double hmin, hmax;
    if (conf==4) { nbin=100; hmin=1.; hmax=4.; } 
    else if (conf==14) { nbin=100; hmin=5.; hmax=10.; }
    else {
      cerr << "[Utilities::TH1FhQ2] Enter valid SBS config!!" << endl;
    }
    TH1F *h = new TH1F(name.c_str(), "Q^{2} Distribution (GeV^{2})", 
		       nbin, hmin, hmax);
    return h;
  }
}

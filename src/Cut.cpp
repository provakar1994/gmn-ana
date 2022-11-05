#include "../include/Cut.h"

namespace cut {

  std::vector<double> hcal_active_area_simu (int nBlk_x=1, int nBlk_y=1) {
    // defines the active area of HCAL using simulation DB
    // indiv. block pos. can be found at SBS-replay/DB_MC/db_sbs.hcal.dat
    std::vector<double> active_area;
    // block positions from DB
    double xHCAL_t_DB_MC = expconst::xHCAL_t_DB_MC;  //m, center of top row blocks (from DB)
    double xHCAL_b_DB_MC = expconst::xHCAL_b_DB_MC;  //m, center of bottom row blocks (from DB)
    double yHCAL_r_DB_MC = expconst::yHCAL_r_DB_MC;  //m, center of right most blocks (from DB)
    double yHCAL_l_DB_MC = expconst::yHCAL_l_DB_MC;  //m, center of left most blocks (from DB)
    // calculate cut limits
    double xHCAL_t = (xHCAL_t_DB_MC + (expconst::hcalblk_h/2.) + expconst::hcalblk_gap_v) + (nBlk_x - 1) * expconst::hcalblk_cTc_v;
    double xHCAL_b = (xHCAL_b_DB_MC - (expconst::hcalblk_h/2.) - expconst::hcalblk_gap_v) - (nBlk_x - 1) * expconst::hcalblk_cTc_v;
    double yHCAL_r = (yHCAL_r_DB_MC + (expconst::hcalblk_w/2.) + expconst::hcalblk_gap_h) + (nBlk_y - 1) * expconst::hcalblk_cTc_h;
    double yHCAL_l = (yHCAL_l_DB_MC - (expconst::hcalblk_w/2.) - expconst::hcalblk_gap_h) - (nBlk_y - 1) * expconst::hcalblk_cTc_h;
    active_area.push_back( xHCAL_t ); 
    active_area.push_back( xHCAL_b );
    active_area.push_back( yHCAL_r );
    active_area.push_back( yHCAL_l );
    return active_area;
  }
  //___________________________________________________________________
  std::vector<double> hcal_active_area_data (int nBlk_x=1, int nBlk_y=1) {
    // defines the active area of HCAL using real data DB
    // indiv. block pos. can be found at SBS-replay/DB/db_sbs.hcal.dat
    std::vector<double> active_area;
    // block positions from DB
    double xHCAL_t_DB = expconst::xHCAL_t_DB;   //m, center of top row blocks (from DB)
    double xHCAL_b_DB = expconst::xHCAL_b_DB;   //m, center of bottom row blocks (from DB)
    double yHCAL_r_DB = expconst::yHCAL_r_DB;   //m, center of right most blocks (from DB)
    double yHCAL_l_DB = expconst::yHCAL_l_DB;   //m, center of left most blocks (from DB)
    // calculate cut limits
    double xHCAL_t = (xHCAL_t_DB + (expconst::hcalblk_h/2.) + expconst::hcalblk_gap_v) + (nBlk_x - 1) * expconst::hcalblk_cTc_v;
    double xHCAL_b = (xHCAL_b_DB - (expconst::hcalblk_h/2.) - expconst::hcalblk_gap_v) - (nBlk_x - 1) * expconst::hcalblk_cTc_v;
    double yHCAL_r = (yHCAL_r_DB + (expconst::hcalblk_w/2.) + expconst::hcalblk_gap_h) + (nBlk_y - 1) * expconst::hcalblk_cTc_h;
    double yHCAL_l = (yHCAL_l_DB - (expconst::hcalblk_w/2.) - expconst::hcalblk_gap_h) - (nBlk_y - 1) * expconst::hcalblk_cTc_h;
    active_area.push_back( xHCAL_t ); 
    active_area.push_back( xHCAL_b );
    active_area.push_back( yHCAL_r );
    active_area.push_back( yHCAL_l );
    return active_area;
  }
  //___________________________________________________________________
  bool inHCAL_activeA (double xHCAL, double yHCAL, std::vector<double> hcal_active_area) {
    // returns "True" if nucleon pos. in HCAL is within active area 
    bool isActive = false;
    // active area dimensions
    double xHCAL_t = hcal_active_area[0];
    double xHCAL_b = hcal_active_area[1];
    double yHCAL_r = hcal_active_area[2];
    double yHCAL_l = hcal_active_area[3];
    isActive = yHCAL>yHCAL_r && yHCAL<yHCAL_l && xHCAL>xHCAL_t && xHCAL<xHCAL_b;
    return isActive;
  } 
  //___________________________________________________________________
  std::vector<double> hcal_safety_margin (double delx_sigma, double dely_sigma, std::vector<double> hcal_active_area) {
    // returns co-ordinates for the safety margin of HCAL using HCAL active area
    std::vector<double> safety_margin;
    double xHCAL_t = hcal_active_area[0] + delx_sigma;  // top margin
    double xHCAL_b = hcal_active_area[1] - delx_sigma;  // bottom margin
    double yHCAL_r = hcal_active_area[2] + dely_sigma;  // right margin
    double yHCAL_l = hcal_active_area[3] - dely_sigma;  // left margin
    safety_margin.push_back( xHCAL_t ); 
    safety_margin.push_back( xHCAL_b );
    safety_margin.push_back( yHCAL_r );
    safety_margin.push_back( yHCAL_l );
    return safety_margin;
  }
  //___________________________________________________________________
  std::vector<double> hcal_safety_margin (double delx_sigma_p, double delx_sigma_n, double dely_sigma, std::vector<double> hcal_active_area) {
    // returns co-ordinates for the safety margin of HCAL using HCAL active area
    std::vector<double> safety_margin;
    double xHCAL_t = hcal_active_area[0] + delx_sigma_p;  // top margin (relevant for proton)
    double xHCAL_b = hcal_active_area[1] - delx_sigma_n;  // bottom margin (relevant for neutron)
    double yHCAL_r = hcal_active_area[2] + dely_sigma;    // right margin
    double yHCAL_l = hcal_active_area[3] - dely_sigma;    // left margin
    safety_margin.push_back( xHCAL_t ); 
    safety_margin.push_back( xHCAL_b );
    safety_margin.push_back( yHCAL_r );
    safety_margin.push_back( yHCAL_l );
    return safety_margin;
  }
  //___________________________________________________________________
  bool inHCAL_fiducial (double xHCAL_exp, double yHCAL_exp, double delx_shift, vector<double> hcal_safety_margin) {
    // returns "True" if expected nucleon pos. in HCAL is within "Fiducial" region
    bool inFidu = false;
    // active area dimensions
    double xHCAL_t = hcal_safety_margin[0];
    double xHCAL_b = hcal_safety_margin[1];
    double yHCAL_r = hcal_safety_margin[2];
    double yHCAL_l = hcal_safety_margin[3];
    // first check whether neutrons are in fiducial region or not
    bool inFidu_n = yHCAL_exp>yHCAL_r && yHCAL_exp<yHCAL_l && xHCAL_exp>xHCAL_t && xHCAL_exp<xHCAL_b;

    // Now, check whether the protons are in fiducial region or not
    // calculate expected xHCAL_exp for proton [considering SBS magnet kick]
    double xHCAL_exp_p = xHCAL_exp - delx_shift;  // "-" sign due to the fact that protons are upbending
    bool inFidu_p = xHCAL_exp_p>xHCAL_t && xHCAL_exp_p<xHCAL_b;  //yHCAL_exp cuts the same for n & p

    // We return "True" if both isFidu_n and isFidu_p are satisfied
    inFidu = inFidu_n && inFidu_p;
    //inFidu = yHCAL_exp>yHCAL_r && yHCAL_exp<yHCAL_l && xHCAL_exp<xHCAL_b && xHCAL_exp_p>xHCAL_t;
    return inFidu;
  } 

}

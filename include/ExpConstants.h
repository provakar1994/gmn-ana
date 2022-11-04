#ifndef EXP_CONSTANTS_H
#define EXP_CONSTANTS_H

#include "TMath.h"

namespace expconst {

  // Detector variables
  // BBCAL
  static const int shcol = 7;
  static const int shrow = 27;
  static const int pscol = 2;
  static const int psrow = 26;
  // HCAL 
  // 1. dimesions found from g4sbs: G4SBSHArmBuilder::MakeHCALV2
  static const int hcalcol = 12;
  static const int hcalrow = 24;
  static const double hcalblk_w = 0.1524;       //m, width of a HCAL block
  static const double hcalblk_h = 0.1524;       //m, height of a HCAL block
  static const double hcalblk_cTc_h = 0.15494;  //m, horizontal center-to-center dist.
  static const double hcalblk_cTc_v = 0.15875;  //m, vertical center-to-center dist.
  static const double hcalblk_gap_h = 0.00254;  //m, horiz. gap bet. two blocks
  static const double hcalblk_gap_v = 0.00635;  //m, vert. gap bet. two blocks
  // 2. x & y pos of blocks sitting at all 4 edges
  static const double xHCAL_t_DB = -2.190625;   //m, center of top row blocks (from SBS-replay/DB)
  static const double xHCAL_b_DB = 1.460625;    //m, center of bottom row blocks (from SBS-replay/DB)
  static const double yHCAL_r_DB = -0.85217;    //m, center of right most blocks (from SBS-replay/DB)
  static const double yHCAL_l_DB = 0.85217;     //m, center of left most blocks (from SBS-replay/DB)
  // --- above => data, below => simu
  static const double xHCAL_t_DB_MC = -2.27563; //m, center of top row blocks (from SBS-replay/DB_MC)
  static const double xHCAL_b_DB_MC = 1.37562;  //m, center of bottom row blocks (from SBS-replay/DB_MC)
  static const double yHCAL_r_DB_MC = -0.85217; //m, center of right most blocks (from SBS-replay/DB_MC)
  static const double yHCAL_l_DB_MC = 0.85217;  //m, center of left most blocks (from SBS-replay/DB_MC)
  /* // 3. Offsets adjusted by looking at deltax and deltay distributions */
  /* static const double hcaloffset_v_data = -0.38; //m, vert. offset of HCAL origin w.r.t DB (data) */
  /* static const double hcaloffset_h_data = 0.15;  //m, horiz. offset of HCAL origin w.r.t DB (data) */
  /* static const double hcaloffset_v_simu = 0.0;   //m, vert. offset of HCAL origin w.r.t DB (simu) */
  /* static const double hcaloffset_h_simu = 0.0;   //m, horiz. offset of HCAL origin w.r.t DB (simu) */
  
  // Constant for the entire experiment
  // target
  static const double tarlen = 15.0;  //cm 
  // LH2
  static const double lh2tarrho = 0.0723;     //g/cc, target density
  static const double lh2cthick = 0.02;       //cm, target cell thickness
  static const double lh2uwallthick = 0.0145; //cm, upstream wall thickness
  static const double lh2dwallthick = 0.015;  //cm, downstream wall thickness
  // LD2
  static const double ld2tarrho = 0.169;      //g/cc, target density

  // magnet
  static const double dipolegap = 1.22;  //m

  // shieldling
  static const double Alrho = 2.7; //g/cc

  // Following quantities vary with configuration
  double ebeam(int config);     //GeV
  double bbtheta(int config);   //deg
  double bbdist(int config);    //m
  double sbstheta(int config);  //deg
  double sbsdist(int config);   //m
  double hcaltheta(int config); //deg
  double hcaldist(int config);  //m
}

// a class for SBS config
class SBSconfig {
 public:

  int    GetSBSconf()       const { return fSBSconf; }
  int    GetSBSmag()        const { return fSBSmag; }
  double GetEbeam()         const { return fEbeam; }
  double GetBBtheta()       const { return fBBtheta; }
  double GetBBtheta_rad()   const { return fBBtheta_rad; }
  double GetBBdist()        const { return fBBdist; }
  double GetSBStheta()      const { return fSBStheta; }
  double GetSBStheta_rad()  const { return fSBStheta_rad; }
  double GetSBSdist()       const { return fSBSdist; }
  double GetHCALtheta()     const { return fHCALtheta; }
  double GetHCALtheta_rad() const { return fHCALtheta_rad; }
  double GetHCALdist()      const { return fHCALdist; }

  // constructor
  SBSconfig(int conf, int sbsmag) {
    fSBSconf       = conf;
    fSBSmag        = sbsmag;
    fEbeam         = expconst::ebeam(conf);
    fBBtheta       = expconst::bbtheta(conf);
    fBBtheta_rad   = expconst::bbtheta(conf)*TMath::DegToRad();
    fBBdist        = expconst::bbdist(conf);
    fSBStheta      = expconst::sbstheta(conf);
    fSBStheta_rad  = expconst::sbstheta(conf)*TMath::DegToRad();
    fSBSdist       = expconst::sbsdist(conf);
    fHCALtheta     = expconst::hcaltheta(conf);
    fHCALtheta_rad = expconst::hcaltheta(conf)*TMath::DegToRad();
    fHCALdist      = expconst::hcaldist(conf);
  }

  // print to screen
  int Print() {
    std::cout << " -------------------------- "                         << std::endl
	      << Form(" SBS Config: %d, "                   , GetSBSconf())   << std::endl
	      << Form(" SBS Magnet Settings: %d (p), "      , GetSBSmag())    << std::endl
    	      << Form(" Beam energy: %0.4f (GeV),"          , GetEbeam())     << std::endl
    	      << Form(" BigBite angle: %0.1f (deg),"        , GetBBtheta())   << std::endl
      	      << Form(" BigBite distance: %0.5f (m),"       , GetBBdist())    << std::endl
    	      << Form(" Super BigBite angle: %0.1f (deg),"  , GetSBStheta())  << std::endl
      	      << Form(" Super BigBite distance: %0.2f (m)," , GetSBSdist())   << std::endl
	      << Form(" HCAL angle: %0.1f (deg),"           , GetHCALtheta()) << std::endl
    	      << Form(" HCAL distance: %0.1f (m)"           , GetHCALdist())  << std::endl
	      << " -------------------------- "                         << std::endl;
    return 0;
  }

 private:
  int    fSBSconf;             // SBS configuration number
  int    fSBSmag;              // SBS magnet settings (%)
  double fEbeam;               // beam energy (better to get this from tree) (GeV)
  double fBBtheta;             // BigBite magnet angle (deg)
  double fBBtheta_rad;         // BigBite magnet angle (rad)
  double fBBdist;              // BigBite magnet distance from target (m)
  double fSBStheta;            // Super BigBite magnet angle (deg)
  double fSBStheta_rad;        // Super BigBite magnet angle (rad)
  double fSBSdist;             // Super BigBite magnet distance from target (m)
  double fHCALtheta;           // HCAL angle (deg)
  double fHCALtheta_rad;       // HCAL angle (rad)
  double fHCALdist;            // HCAL distance from target (m)
};

#endif

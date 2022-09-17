#ifndef UTIL_SBS_CONFIG_H
#define UTIL_SBS_CONFIG_H

#include "../src/ExpConstants.cpp"

// a class for SBS config
class SBSconfig {
 public:

  int    GetSBSconf()      const { return fSBSconf; }
  int    GetSBSmag()       const { return fSBSmag; }
  double GetEbeam()        const { return fEbeam; }
  double GetBBtheta()      const { return fBBtheta; }
  double GetBBtheta_rad()  const { return fBBtheta_rad; }
  double GetBBdist()       const { return fBBdist; }
  double GetSBStheta()     const { return fSBStheta; }
  double GetSBStheta_rad() const { return fSBStheta_rad; }
  double GetSBSdist()      const { return fSBSdist; }
  double GetHCALdist()     const { return fHCALdist; }

  // constructor
  SBSconfig(int conf, int sbsmag) {
    fSBSconf      = conf;
    fSBSmag       = sbsmag;
    fEbeam        = expconst::ebeam(conf);
    fBBtheta      = expconst::bbtheta(conf);
    fBBtheta_rad  = expconst::bbtheta(conf)*TMath::DegToRad();
    fBBdist       = expconst::bbdist(conf);
    fSBStheta     = expconst::sbstheta(conf);
    fSBStheta_rad = expconst::sbstheta(conf)*TMath::DegToRad();
    fSBSdist      = expconst::sbsdist(conf);
    fHCALdist     = expconst::hcaldist(conf);
  }

  // print to screen
  int Print() {
    std::cout << " -------------------------- "                         << std::endl
	      << Form(" SBS Config: %d, "                   , GetSBSconf())  << std::endl
	      << Form(" SBS Magnet Settings: %d (p), "      , GetSBSmag())   << std::endl
    	      << Form(" Beam energy: %0.4f (GeV),"          , GetEbeam())    << std::endl
    	      << Form(" BigBite angle: %0.1f (deg),"        , GetBBtheta())  << std::endl
      	      << Form(" BigBite distance: %0.5f (m),"       , GetBBdist())   << std::endl
    	      << Form(" Super BigBite angle: %0.1f (deg),"  , GetSBStheta()) << std::endl
      	      << Form(" Super BigBite distance: %0.2f (m)," , GetSBSdist())  << std::endl
    	      << Form(" HCAL distance: %0.1f (m)"           , GetHCALdist()) << std::endl
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
  double fHCALdist;            // HCAL distance from target (m)
};

#endif

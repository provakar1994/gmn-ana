#ifndef UTIL_SBS_CONFIG_H
#define UTIL_SBS_CONFIG_H

// a struct for SBS config

typedef struct SBSconfig {
  int    sbsconf;             // SBS configuration number
  double Ebeam;               // beam energy (better to get this from tree)
  double BBtheta;             // BigBite magnet angle
  double BBdist;              // BigBite magnet distance from target
  double SBStheta;            // Super BigBite magnet angle
  double SBSdist;             // Super BigBite magnet distance from target
  double HCALdist;            // HCAL distance from target

  // constructor
  SBSconfig():
    Ebeam(0), BBtheta(0), BBdist(0), SBStheta(0), SBSdist(0), HCALdist(0)
  {}

  int    GetSBSconf()  const { return sbsconf; }
  double GetEbeam()    const { return Ebeam; }
  double GetBBtheta()  const { return BBtheta; }
  double GetBBdist()   const { return BBdist; }
  double GetSBStheta() const { return SBStheta; }
  double GetSBSdist()  const { return SBSdist; }
  double GetHCALdist() const { return HCALdist; }

  // print to screen
  int Print() {
    std::cout << " -------------------------- "                         << std::endl
	      << Form(" SBS Config: %d, "                   , sbsconf)  << std::endl
    	      << Form(" Beam energy: %0.4f (GeV),"          , Ebeam)    << std::endl
    	      << Form(" BigBite angle: %0.1f (deg),"        , BBtheta)  << std::endl
      	      << Form(" BigBite distance: %0.5f (m),"       , BBdist)   << std::endl
    	      << Form(" Super BigBite angle: %0.1f (deg),"  , SBStheta) << std::endl
      	      << Form(" Super BigBite distance: %0.2f (m)," , SBSdist)  << std::endl
    	      << Form(" HCAL distance: %0.1f (m)"           , HCALdist) << std::endl
	      << " -------------------------- "                         << std::endl;
    return 0;
  }
} SBSconfig_t;

#endif

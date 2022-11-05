#ifndef CUT_H
#define CUT_H

#include <iostream>
#include "TString.h"
#include "ExpConstants.h"

namespace cut {


  
  /* ##########################
     ##  HCAL Geometry Cuts  ##
     ########################## */

  // Defines active area of HCAL excluding desired no. of blocks from all 4 sides
  // Uses ** Simulation ** DB to calculate the active region
  // Default: Excludes 1 block from all 4 sides of HCAL
  std::vector<double> hcal_active_area_simu (int nBlk_x ,  // No. of blocks to exclude form top and bottom (Default=1) 
					     int nBlk_y);  // No. of blocks to exclude form left and right (Default=1)

  // Defines active area of HCAL excluding desired no. of blocks from all 4 sides
  // Uses ** Real Data ** DB to calculate the active region
  // Default: Excludes 1 block from all 4 sides of HCAL
  std::vector<double> hcal_active_area_data (int nBlk_x,   // No. of blocks to exclude form top and bottom (Default=1)
					     int nBlk_y);  // No. of blocks to exclude form left and right (Default=1)

  // Returns "True" if nucleon pos. in HCAL is within active area
  bool isHCAL_activeA (double xHCAL,                           // vertical (x) pos of recoil N at HCAL
		       double yHCAL,                           // horizontal (y) pos of recoil N at HCAL
		       std::vector<double> hcal_active_area);  // HCAL active area co-ordinates

  // Defines safety region within HCAL for expected nucleon positions
  std::vector<double> hcal_safety_margin (double delx_sigma,                      // vertical (x) pos of recoil N at HCAL
					  double dely_sigma,                      // horizontal (y) pos of recoil N at HCAL
					  std::vector<double> hcal_active_area);  // HCAL active area co-ordinates
  // Defines safety region within HCAL for expected nucleon positions (Overloading)
  std::vector<double> hcal_safety_margin (double delx_sigma_p,                    // vertical (x) pos of recoil p at HCAL
					  double delx_sigma_n,                    // vertical (x) pos of recoil n at HCAL
					  double dely_sigma,                      // horizontal (y) pos of recoil N at HCAL
					  std::vector<double> hcal_active_area);  // HCAL active area co-ordinates

  // Returns "True" if expected nucleon pos. in HCAL is within "Fiducial" region
  bool inHCAL_fiducial (double xHCAL_exp,                        // m, expected vert. (x) pos of recoil N at HCAL
			double yHCAL_exp,                        // m, expected horiz. (y) pos of recoil N at HCAL
			double delx_shift,                       // m, amount of SBS magnet kick
			std::vector<double> hcal_safety_margin); // HCAL safety margin co-ordinates

  /* ######################## */
}

#endif

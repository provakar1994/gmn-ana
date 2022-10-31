/* 
   gmn-ana.h will include all the gmn-ana libraries in the most
   efficient way. Just including this file in an analysis script
   will be enough to get access to all the gmn-ana libraries.
   -----
   P. Datta <pdbforce@jlab.org> Created 09-17-2022
*/

#include "Constants.h"                // namespace constant
#include "../src/Cut.cpp"             // namespace cut
#include "../src/Utilities.cpp"       // namespace util_pd
#include "../src/ExpConstants.cpp"    // namespace expconst & class SBSconfig
#include "../src/SetROOTVar.cpp"      // namespace setrootvar
#include "../src/KinematicVar.cpp"    // namespace kine

/* --- List of gmn-ana libraries --- */
// Cut.h          : namespace Cut
// Constants.h    : namespace constant
// KinematicVar.h : namespace kine
// ExpConstants.h : namespace expconst & class SBSconfig
// SetROOTVar.h   : namespace setrootvar
// Utilities.h    : namespace util_pd

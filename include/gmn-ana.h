/* 
   gmn-ana.h will include all the gmn-ana libraries in the most
   efficient way. Just including this file in an analysis script
   will be enough to get access to all the gmn-ana libraries.
   -----
   P. Datta <pdbforce@jlab.org> Created 09-17-2022
*/

#include "Constants.h"                // namespace constant
#include "SBSconfig.h"                // class sbsconf
#include "../src/SetROOTVar.cpp"      // namespace setrootvar
#include "../src/KinematicVar.cpp"    // namespace kine

/* --- List of gmn-ana libraries --- */
// SBSconfig.h    : class SBSconfig
// Constants.h    : namespace constant
// KinematicVar.h : namespace kine
// ExpConstants.h : namespace expconst
// SetROOTVar.h   : namespace setrootvar

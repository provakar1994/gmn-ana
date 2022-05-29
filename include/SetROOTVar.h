#ifndef SET_ROOT_VAR_H
#define SET_ROOT_VAR_H

#include "TTree.h"

namespace setvar {

  // PS cluster
  /* void setpsclusvar(); */
  /* void setpsclustimevar(TTree *T, double ps_atimeblk); */

  /* void setshclvar(); */

  /* void sethcalclvar(); */

  void setbranch(TTree *T,
		 std::string prefix,
		 std::vector<std::string> suffix,
		 std::vector<double*> memory);

}

#endif

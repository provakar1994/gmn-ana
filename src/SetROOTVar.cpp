#include "../include/SetROOTVar.h"

namespace setvar {

  // PS Cluster

  // void setpsclustimevar( TTree *T, double &ps_atimeblk){
  //   T->SetBranchStatus("bb.ps.atimeblk",1);
  //   T->SetBranchAddress("bb.ps.atimeblk",ps_atimeblk);
  // }

  void setbranch(TTree *T,
		 std::string prefix,
		 std::vector<std::string> suffix,
		 std::vector<double*> memory)
  {
    std::string branchname;
    for(int i=0;i<suffix.size();i++){
      branchname = prefix + std::string(".") + suffix[i];
      T->SetBranchStatus(branchname.c_str(),1);
      T->SetBranchAddress(branchname.c_str(),memory[i]);
    }
  }

}

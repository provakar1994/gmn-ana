#include "../include/SetROOTVar.h"

namespace setrootvar {

  void setbranch(TTree *T,
		 std::string prefix,
		 std::string suffix,
		 void* memory)
  {
    std::string branchname;
    branchname = prefix + std::string(".") + suffix;
    T->SetBranchStatus(branchname.c_str(),1);
    T->SetBranchAddress(branchname.c_str(),memory);
  }
  //--------------------------------------------
  void setbranch(TTree *T,
		 std::string prefix,
		 std::vector<std::string> suffix,
		 std::vector<void*> memory)
  {
    if(suffix.size()!=memory.size()){
      std::cerr << "Contaier size of branch suffix != Container size of branch memory" std::endl;
      throw;
    }
    std::string branchname;
    for(int i=0;i<suffix.size();i++){
      branchname = prefix + std::string(".") + suffix[i];
      T->SetBranchStatus(branchname.c_str(),1);
      T->SetBranchAddress(branchname.c_str(),memory[i]);
    }
  }
  //--------------------------------------------
  void setbranch(TTree *T,
  		 std::string prefix,
  		 std::vector<std::string> suffix,
  		 std::vector<void*> memory,
  		 int ndatavarpos)
  {
    if(suffix.size()!=memory.size()){
      std::cerr << "Contaier size of branch suffix != Container size of branch memory" std::endl;
      throw;
    }
    std::string branchname;
    for(int i=0;i<suffix.size();i++){
      if(i!=ndatavarpos){
  	branchname = prefix + std::string(".") + suffix[i];
  	T->SetBranchStatus(branchname.c_str(),1);
  	T->SetBranchAddress(branchname.c_str(),memory[i]);
      }else{
  	branchname = std::string("Ndata.")+prefix+std::string(".")+suffix[i];
  	T->SetBranchStatus(branchname.c_str(),1);
  	T->SetBranchAddress(branchname.c_str(),memory[i]);
      }
    }
  }
  //--------------------------------------------
  void setbranch(TTree *T,
  		 std::string prefix,
  		 std::vector<std::string> suffix,
  		 std::vector<void*> memory,
  		 std::vector<int> ndatapos)
  {
    if(suffix.size()!=memory.size()){
      std::cerr << "Contaier size of branch suffix != Container size of branch memory" std::endl;
      throw;
    }
    std::string branchname;
    for(int i=0;i<suffix.size();i++){
      if(std::find(ndatapos.begin(),ndatapos.end(),i)!=ndatapos.end()){
  	branchname = std::string("Ndata.")+prefix+std::string(".")+suffix[i];
  	T->SetBranchStatus(branchname.c_str(),1);
  	T->SetBranchAddress(branchname.c_str(),memory[i]);
      }else{
 	branchname = prefix + std::string(".") + suffix[i];
  	T->SetBranchStatus(branchname.c_str(),1);
  	T->SetBranchAddress(branchname.c_str(),memory[i]);
      }
    }
  }
  //--------------------------------------------
}

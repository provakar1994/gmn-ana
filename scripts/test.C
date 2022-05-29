#include <iostream>
#include <vector>

#include "../src/ExpConstants.cpp"
#include "../src/KinematicVar.cpp"
//#include "../src/SetROOTVar.cpp"

int main()
{
  //std::cout << kine::W(2,2,kine::Q2(2,2,0.2)) << std::endl;

  //TTree *T;
  
  //bbsh clus var
  std::vector<std::string> shclvar = {"e","x","y","rowblk","colblk"};
  double eSH, xSH, ySH, rblkSH, cblkSH;
  vector<double*> shclvar_mem = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH};

  std::string added = "bb.sh" + std::string(".") + shclvar[0];
  //std::cout << added.c_str() << std::endl;

  //setvar::setbranch(T,"bb.sh",shclvar,shclvar_mem);

  // double x=0., y[2]={1.,2}, z=2.;
  // std::string xstr[3]={"a","b","gh"};
  // vector<double*> xa = {&x,&z};
  //xa.push_back(&x); xa.push_back(&z);
  //void* ya=&y;
  //void* za[2]={&x,&z};

  //void* abc = xa;

  std::cout << shclvar.size() << std::endl;

  // std::vector<int> vint(12);

  // std::vector<char> vstr('a');
  
  // std::vector<void*> vec;
  // vec.push_back(&y);

  //std::cout << shclvar[3].c_str() << " " << shclvar_mem[3] << std::endl;
  
  return 0;
}

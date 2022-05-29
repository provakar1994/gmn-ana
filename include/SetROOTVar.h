#ifndef SET_ROOT_VAR_H
#define SET_ROOT_VAR_H

#include "TTree.h"

namespace setrootvar {

  void setbranch(TTree *T, //use to settreebranch for one variable
		 std::string prefix,
		 std::string suffix,
		 void* memory);
  void setbranch(TTree *T, //use to settreebranch for multiple variables
		 std::string prefix,
		 std::vector<std::string> suffix,
		 std::vector<void*> memory);
  void setbranch(TTree *T, //use if one of the variables is Ndata
		 std::string prefix,
		 std::vector<std::string> suffix,
		 std::vector<void*> memory,
		 int ndatavarpos);
  void setbranch(TTree *T, //use if there are multiple Ndata variables
		 std::string prefix,
		 std::vector<std::string> suffix,
		 std::vector<void*> memory,
		 std::vector<int> ndatapos);

}

#endif

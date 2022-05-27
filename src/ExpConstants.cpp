#include "../include/ExpKineConstants.h"

namespace expconst {

  double ebeam(int config){
    if(config==1)
      return 1.916;
    else if(config==4)
      return 3.7278;
    else if(config==7)
      return 7.906;
    else if(config==11)
      return 9.91;
    else if(config==14 || config==8)
      return 5.965;
    else if(config==9)
      return 4.013;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

  double bbtheta(int config){
    if(config==1)
      return 51.0;
    else if(config==4)
      return 36.0;
    else if(config==7)
      return 40.0;
    else if(config==11)
      return 42.0;
    else if(config==14)
      return 46.5;
    else if(config==8)
      return 26.5;
    else if(config==9)
      return 49.0;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

  double bbdist(int config){
    if(config==1)
      return 1.8518;
    else if(config==4)
      return 1.7988;
    else if(config==7)
      return 1.84896;
    else if(config==11)
      return 1.55146;
    else if(config==14)
      return 1.84787;
    else if(config==8)
      return 1.97473;
    else if(config==9)
      return 1.550;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

  double sbstheta(int config){
    if(config==1)
      return 33.5;
    else if(config==4)
      return 31.9;
    else if(config==7)
      return 16.1;
    else if(config==11)
      return 13.3;
    else if(config==14)
      return 17.3;
    else if(config==8)
      return 29.9;
    else if(config==9)
      return 22.5;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

  double sbsdist(int config){
    if(config==1||config==4||config==7||config==11
       ||config==14||config==8||config==9)
      return 2.25;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

  double hcaldist(int config){
    if(config==1)
      return 13.5;
    else if(config==4||config==8||config==9)
      return 11.0;
    else if(config==7||config==14)
      return 14.0;
    else if(config==11)
      return 14.5;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

} //::expconst



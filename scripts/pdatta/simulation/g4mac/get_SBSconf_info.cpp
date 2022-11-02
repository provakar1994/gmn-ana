#include <iostream>

#include "../../../../include/gmn-ana.h"

void get_SBSconf_info (int conf=4, int sbsmag=30) {

  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

}

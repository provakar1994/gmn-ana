/*
  In simulation HCAL coordinate system has the following convention:
*/

#include <iostream>

const int kNcols = 12; // # of HCAL columns
const int kNrows = 24; // # of HCAL rows
const double hgap = 0.15494; //m, horizontal gap between the centers of two modules
const double vgap = 0.15875; //m, vertical gap between the centers of two modules

int main()
{
  double voffset = 0.45; //m
  double hoffset = 0.0; //m


  /* Co-ordinate system convention (same as real data)
  1. +x points to hall floor
  2. +y points towards beamline (HCAL is on beam right)
  3. +z points along particle trajectory
  */
  // Module # 0 of HCAL is situated at the top right corner while being
  // viewed from front. Let's calculate its co-ordinates.
  double y_topRB = -( voffset + (((double)kNrows-1.)/2.)*vgap );
  double x_topRB = -( hoffset + (((double)kNcols-1.)/2.)*hgap );

  //generating x-positions (vertical positions)
  std::cout << "sbs.hcal.xpos = " << std::endl;
  for(int row=0; row<kNrows; row++){
    for(int col=0; col<kNcols; col++){
      std::cout << y_topRB+row*vgap << " ";
    }
    std::cout << std::endl;
  }
  
  std::cout << std::endl << std::endl;

  //generating y-positions (horizontal positions)
  std::cout << "sbs.hcal.ypos = " << std::endl;
  for(int row=0; row<kNrows; row++){
    for(int col=0; col<kNcols; col++){
      std::cout << x_topRB+col*hgap << " ";
    }
    std::cout << std::endl;
  }

  /* Co-ordinate system convention
  1. +y points to hall roof
  2. +x points towards beamline (HCAL is on beam right)
  3. +z points along particle trajectory
  */

  // Module # 0 of HCAL is situated at the top right corner while being
  // viewed from front. Let's calculate its co-ordinates.
  // double y_topRB = voffset + (((double)kNrows-1.)/2.)*vgap;
  // double x_topRB = hoffset - (((double)kNcols-1.)/2.)*hgap;

  // //generating y-positions (vertical positions)
  // std::cout << "sbs.hcal.ypos = " << std::endl;
  // for(int row=0; row<kNrows; row++){
  //   for(int col=0; col<kNcols; col++){
  //     std::cout << y_topRB-row*vgap << " ";
  //   }
  //   std::cout << std::endl;
  // }
  
  // std::cout << std::endl << std::endl;

  // //generating x-positions (horizontal positions)
  // std::cout << "sbs.hcal.xpos = " << std::endl;
  // for(int row=0; row<kNrows; row++){
  //   for(int col=0; col<kNcols; col++){
  //     std::cout << x_topRB+col*hgap << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // Module # 0 of HCAL is situated at the top left corner while being
  // viewed from front. Let's calculate its co-ordinates.
  // double y_topRB = voffset + (((double)kNrows-1.)/2.)*vgap;
  // double x_topRB = hoffset + (((double)kNcols-1.)/2.)*hgap;

  // //generating y-positions (vertical positions)
  // std::cout << "sbs.hcal.ypos = " << std::endl;
  // for(int row=0; row<kNrows; row++){
  //   for(int col=0; col<kNcols; col++){
  //     std::cout << y_topRB-row*vgap << " ";
  //   }
  //   std::cout << std::endl;
  // }
  
  // std::cout << std::endl << std::endl;

  // //generating x-positions (horizontal positions)
  // std::cout << "sbs.hcal.xpos = " << std::endl;
  // for(int row=0; row<kNrows; row++){
  //   for(int col=0; col<kNcols; col++){
  //     std::cout << x_topRB-col*hgap << " ";
  //   }
  //   std::cout << std::endl;
  // }


  return 0;
}

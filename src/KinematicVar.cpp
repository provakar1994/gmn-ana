#include "../include/KinematicVar.h"

namespace kine {

  //--------------------------------------------
  double pcentral(double ebeam, double etheta, std::string Ntype) {
    double temp = 0.;
    if (Ntype.compare("p")) 
      temp = ebeam/(1. + (ebeam/constant::Mp)*(1.0 - cos(etheta)));
    else if (Ntype.compare("n")) 
      temp = ebeam/(1. + (ebeam/constant::Mn)*(1.0 - cos(etheta)));
    else
      std::cerr << "[KinematicVar::pcentral] Enter a valid nucleon type! **!**" << std::endl;
    return temp;
  }
  //--------------------------------------------
  double etheta(TLorentzVector Peprime) {
    return acos(Peprime.Pz() / Peprime.E());
  }
  //--------------------------------------------
  double ephi(TLorentzVector Peprime) {
    return atan2(Peprime.Py(), Peprime.Px());
  }
  //--------------------------------------------
  void SetPN(std::string Ntype, TLorentzVector &PN) {
    if (Ntype.compare("p")) 
      PN.SetPxPyPzE(0., 0., 0., constant::Mp);
    else if (Ntype.compare("n")) 
      PN.SetPxPyPzE(0., 0., 0., constant::Mn);
    else
      std::cerr << "[KinematicVar::pcentral] Enter a valid nucleon type! **!**" << std::endl;
  } 
  //--------------------------------------------
  double pN_expect(double nu, std::string Ntype) {
    if (Ntype.compare("p"))                      
      return sqrt(pow(nu, 2.) + 2. * constant::Mp * nu);
    else if (Ntype.compare("n"))      
      return sqrt(pow(nu, 2.) + 2. * constant::Mn * nu);
    else {
      std::cerr << "[KinematicVar::pN_expect] Enter a valid nucleon type! **!**" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  TVector3 qVect_unit(double Ntheta, double Nphi) {
    TVector3 pNhat(sin(Ntheta) * cos(Nphi), sin(Ntheta) * sin(Nphi), cos(Ntheta));
    return pNhat;
  }
  //--------------------------------------------
  void SetHCALaxes(double sbstheta_rad,                           // SBS angle in radian 
		   vector<TVector3> &HCAL_axes) {
    TVector3 HCAL_zaxis(sin(-sbstheta_rad),0,cos(-sbstheta_rad)); // Clock-wise rotation about Y axis
    TVector3 HCAL_xaxis(0,-1,0);                                  // -Y axis of Hall CoS = X axis of HCAL CoS
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
    HCAL_axes.push_back(HCAL_xaxis);
    HCAL_axes.push_back(HCAL_yaxis);
    HCAL_axes.push_back(HCAL_zaxis);
  }
  //--------------------------------------------
  void GetxyHCALexpect(TVector3 vertex, TVector3 pNhat, TVector3 HCAL_origin, 
		       vector<TVector3> HCAL_axes, vector<double> &xyHCALexpect) {
    /* This function calculates the expected vertical (x) and horizontal (y) positions
     of the recoil nucleon at the face of HCAL. 
     input:
     1. vertex       : vertex vector [in Hall CoS?? It must be but haven't confirmed], 
     2. pNhat        : projected q vector, 
     3. HCAL_origin  : HCAL origin vector [in Hall CoS], 
     4. HCAL_axes    : HCAL CoS axes [in Hall CoS]
     output:
     1. xyHCALexpect : expected x and y positions
    */
    // Intersection of a ray with a plane
    double sintersect = (HCAL_origin - vertex).Dot(HCAL_axes[2]) / (pNhat.Dot(HCAL_axes[2]));
    // ray from Hall origin onto the face of HCAL where the nucleon hit
    TVector3 HCAL_intersect = vertex + sintersect*pNhat; 

    double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot(HCAL_axes[0]);
    double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot(HCAL_axes[1]);

    xyHCALexpect.push_back(xexpect_HCAL);
    xyHCALexpect.push_back(yexpect_HCAL);
  }
  //--------------------------------------------
  double Q2(double ebeam, double eeprime, double etheta) {
    return 2.0*ebeam*eeprime*(1.0-cos(etheta));
  }
  //--------------------------------------------
  double W2(double ebeam, double eeprime, double Q2, std::string Ntype) {
    double temp = 0.;
    if (Ntype.compare("p")) 
      temp = pow(constant::Mp,2.0) + 2.0*constant::Mp*(ebeam-eeprime) - Q2;
    else if (Ntype.compare("n")) 
      temp = pow(constant::Mn,2.0) + 2.0*constant::Mn*(ebeam-eeprime) - Q2;
    else
      std::cerr << "[KinematicVar::W2] Enter a valid nucleon type! **!**" << std::endl;
    return temp;
  }
  //--------------------------------------------
  double W(double ebeam, double eeprime, double Q2, std::string Ntype) {
    return max(0., sqrt(kine::W2(ebeam, eeprime, Q2, Ntype)));
  }
  //--------------------------------------------
  double Luminosity(double ibeam, std::string targetType) {
    double lumi = 0.;
    if (targetType.compare("LH2"))
      lumi = ((ibeam/constant::qe)*expconst::tarlen*expconst::lh2tarrho*(constant::N_A/constant::D2_Amass));
    else if (targetType.compare("LD2"))
      lumi = ((ibeam/constant::qe)*expconst::tarlen*expconst::ld2tarrho*(constant::N_A/constant::D2_Amass));
    else
      std::cerr << "[KinematicVar::Luminosity] Enter a valid target type! **!**" << std::endl;
    return lumi;
  }

} //::kine
